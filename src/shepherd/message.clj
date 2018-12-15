(ns shepherd.message
  (:require
   [clojure.edn :as edn]
   [clojure.walk :as walk]
   [clojure.java.io :as io]
   [manifold.bus :as bus]
   [manifold.stream :as stream]
   [manifold.deferred :as defer]
   [aleph.http :as http]
   [cheshire.core :as json]
   [cheshire.factory :as factory]
   [kinsky.client :as kafka]
   [polaris.core :as polaris]
   [taoensso.timbre :as log]
   [ring.middleware.params :as params]
   [ring.middleware.resource :as resource]
   [ring.middleware.keyword-params :as keyword])
  (:import
    [iff ChunkReader ChunkWriter]
    org.apache.kafka.common.serialization.ByteArrayDeserializer
    org.apache.kafka.common.serialization.ByteArraySerializer))

(def poll-interval Long/MAX_VALUE)
(def non-numeric-factory
  (factory/make-json-factory {:allow-non-numeric-numbers true}))

(defn resource
  [path]
  (-> path
      io/resource
      io/input-stream
      io/reader))

(defn read-config
  [path]
  (edn/read
   (java.io.PushbackReader.
    (resource path))))

(defn producer-config
  [config]
  {:bootstrap.servers (:host config)})


(defn byte-array-serializer
  "Kafka's own byte[] serializer."
  []
  (ByteArraySerializer.))

(defn byte-array-deserializer
  "Kafka's own byte[] deserializer"
  []
  (ByteArrayDeserializer.))


(defn write-chunk
  [^String type ^"[B" body ^java.io.InputStream stream]
  (.write (ChunkWriter. type body) stream false))

(defn serialize-to-chunks
  "Serialize an Agent message dict with optional :blobs list to a byte[] of chunks."
  [dict]
  (let [blobs (get dict :blobs [])
        agent-msg (dissoc dict :blobs)
        stream (java.io.ByteArrayOutputStream.)]
    (write-chunk "JSON" (.getBytes (json/generate-string agent-msg)) stream)
    (doseq [blob blobs] (write-chunk "BLOB" blob stream))
    (.toByteArray stream)))

(defn read-chunks
  [^java.io.InputStream stream]
  (loop [chunks (ChunkReader/readAll stream false)
         agent-msg {}
         blobs []]
    (let [chunk (first chunks)]
      (if chunk
        (case (.-chunkType chunk)
          "JSON"
          (recur
            (rest chunks)
            (if (empty? agent-msg)
              (json/parse-string (String. (.-body chunk) "UTF-8") true)
              agent-msg)
            blobs)
          "BLOB"
          (recur
            (rest chunks)
            agent-msg
            (conj blobs (.-body chunk)))
          )
        (assoc agent-msg :blobs blobs)))))

(defn deserialize-from-chunks
  "Deserialize an Agent message map with :blobs list from a byte[] of chunks."
  [bytes]
  (read-chunks (java.io.ByteArrayInputStream. bytes)))


(defn boot-producer
  [config]
  (kafka/producer
   (producer-config config)
   (kafka/keyword-serializer)
   (byte-array-serializer)))

(defn send!
  [producer topic message]
  (kafka/send!
   producer
   {:topic topic
    :value (serialize-to-chunks message)}))

(defn consumer-config
  [config]
  {:bootstrap.servers (:host config)
   :enable.auto.commit "true"
   :auto.commit.interval.ms "1000"
   :group.id (get config :group-id "flow")
   :auto.offset.reset "latest"})

(defn boot-consumer
  [config]
  (kafka/consumer
   (consumer-config config)
   (kafka/keyword-deserializer)
   (byte-array-deserializer)))

(defn handle-message
  [state bus producer handle record]
  (try
    (if (= (first record) :by-topic)
      (let [topics (last record)]
        (doseq [[topic messages] topics]
          (doseq [message messages]
            (let [payload (:value message)
                  agent-msg (deserialize-from-chunks payload)
                  msg {topic (dissoc agent-msg :blobs)}]
              (log/info (get-in msg ["shepherd-receive" :event] "") msg)
              (handle bus producer topic agent-msg)
              (swap! state assoc :last-message msg)
              (bus/publish! bus topic (json/generate-string msg)))))))
    (catch Exception e
      (log/error (.getMessage e))
      (.printStackTrace e))))

(defn poll!
  [consumer]
  (try
    (kafka/poll! consumer poll-interval)
    (catch Exception e
      (log/error (.getMessage e)))))

(defn consume
  [consumer handle]
  (binding [factory/*json-factory* non-numeric-factory]
    (loop [records (poll! consumer)]
      (when-not (empty? records)
        (doseq [record records]
          (handle record)))
      (recur (poll! consumer)))))

(defn boot-kafka
  [state config]
  (let [bus (bus/event-bus)
        producer (boot-producer config)
        consumer (boot-consumer config)
        handle (get config :handle-message (fn [_ _]))]
    (doseq [topic (:subscribe config)]
      (log/info "subscribing to topic" topic)
      (kafka/subscribe! consumer topic))
    {:bus bus
     :producer producer
     :consumer
     (defer/future
       (consume
        consumer
        (partial handle-message state bus producer handle)))}))

(defn default-handle-client
  [state conn message]
  (condp = (:event message)
    "INITIALIZE"
    (let [last-message (get @(:state state) :last-message)]
      (stream/put! conn (json/generate-string last-message)))

    (let [producer (:producer state)]
      (send! producer "flow" message))))

(defn boot
  [{:keys [kafka] :as config}]
  (let [state (atom {:last-message {}})
        {:keys [bus producer consumer]} (boot-kafka state kafka)
        handle-client (get config :handle-client default-handle-client)]
    {:config config
     :bus bus
     :producer producer
     :consumer consumer
     :state state
     :handle-client handle-client}))

(defn index-handler
  [state]
  (fn [request]
    {:status 200
     :headers {"content-type" "text/html"}
     :body (slurp "resources/public/index.html")}))

(defn connect-websocket
  [request]
  (defer/catch
      (http/websocket-connection request)
      (fn [_])))

(def non-websocket
  {:status 400
   :headers {"content-type" "application/text"}
   :body "must connect using websocket request"})

(defn websocket-handler
  [{:keys [bus config handle-client] :as state}]
  (fn [request]
    (defer/let-flow [conn (connect-websocket request)]
      (if-not conn
        non-websocket
        (do
          (stream/connect
           (bus/subscribe
            bus
            (get-in config [:kafka :event-topic]))
           conn)
          (stream/consume
           (fn [from-client]
             (let [from (json/parse-string from-client true)]
               (log/info "client :" from)
               (handle-client state conn from)))
           conn))))))

(defn make-flow-routes
  [state]
  [["/" :index (#'index-handler state)]
   ["/ws" :websocket (#'websocket-handler state)]])

(defn start
  [state]
  (let [config (:config state)
        make-routes (get config :routes (fn [_] []))
        flow-routes (make-flow-routes state)
        other-routes (make-routes state)
        routes (polaris/build-routes (concat flow-routes other-routes))
        router (polaris/router routes)
        app (-> router
                (resource/wrap-resource "public")
                (keyword/wrap-keyword-params)
                (params/wrap-params))]
    (http/start-server app {:port (:port config)})))
