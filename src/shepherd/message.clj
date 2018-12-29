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
    [chunk ChunkReader ChunkWriter]
    org.apache.kafka.common.serialization.ByteArrayDeserializer
    org.apache.kafka.common.serialization.ByteArraySerializer))

(set! *warn-on-reflection* true) ; DEBUG

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


(defn write-chunk!
  "Write a chunk to a stream."
  [^String type ^"[B" body ^java.io.OutputStream stream]
  (.write (ChunkWriter. type body) stream false))

(defn serialize-to-chunks
  "Serialize an Agent message map with an optional :blobs list to a Kafka
  message payload in chunk format."
  [agent-message]
  (let [blobs (get agent-message :blobs [])
        agent-msg (dissoc agent-message :blobs)
        stream (java.io.ByteArrayOutputStream.)]
    (write-chunk! "JSON" (.getBytes (json/generate-string agent-msg)) stream)
    (doseq [blob blobs]
      (write-chunk! "BLOB" blob stream))
    (.toByteArray stream)))

(defn message-with-optional-blobs
  [agent-msg blobs]
  (if (empty? blobs)
    agent-msg
    (assoc agent-msg :blobs blobs)))

(defn read-json-chunk
  [^ChunkWriter chunk]
  (let [body (String. (.-body chunk) "UTF-8")]
    (json/parse-string body true)))

(defn decode-json-chunk
  [previous-msg ^ChunkWriter chunk]
  (if (empty? previous-msg) (read-json-chunk chunk) previous-msg))

(defn read-chunks!
  "Read an Agent message map with an optional :blobs list from a stream of chunks."
  [^java.io.InputStream payload-stream]
  (loop [chunks (ChunkReader/readAll payload-stream false)
         agent-msg {}
         blobs []]
    (let [^ChunkWriter chunk (first chunks)
          more (rest chunks)]
      (case (if chunk (.-chunkType chunk) :done)
        :done
        (message-with-optional-blobs agent-msg blobs)
        "JSON"
        (recur more (decode-json-chunk agent-msg chunk) blobs)
        "BLOB"
        (recur more agent-msg (conj blobs (.-body chunk)))
        (recur more agent-msg blobs)))))

(defn deserialize-from-chunks
  "Deserialize an Agent message map with an optional :blobs list from a Kafka
  message payload in chunk format."
  [^"[B" payload-bytes]
  (read-chunks! (java.io.ByteArrayInputStream. payload-bytes)))


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
                  agent-message (deserialize-from-chunks payload)
                  num-blobs (count (:blobs agent-message))
                  blob-note (if (pos? num-blobs) (str "+ " num-blobs " BLOBs") "")
                  agent-msg (dissoc agent-message :blobs)
                  topic-msg-pair {topic agent-msg}]
              (log/info (:event agent-msg "") topic-msg-pair blob-note)
              (handle bus producer topic agent-message)
              (swap! state assoc :last-message topic-msg-pair)
              (bus/publish! bus topic (json/generate-string topic-msg-pair)))))))
    (catch Exception e
      (log/error e))))

(defn poll!
  [consumer]
  (try
    (kafka/poll! consumer poll-interval)
    (catch Exception e
      (log/error e))))

(defn consume
  [consumer handle]
  (binding [factory/*json-factory* non-numeric-factory]
    ; TODO(jerry): Revisit unpacking the poll! result.
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
