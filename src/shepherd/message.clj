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
  "Load a resource at the given path."
  [path]
  (-> path
      io/resource
      io/input-stream
      io/reader))

(defn read-config
  "Read a configuration map from the given path."
  [path]
  (edn/read
   (java.io.PushbackReader.
    (resource path))))

(defn producer-config
  "Generate a config for a kafka producer, given a map of config values.
     config - map of configuration values for the kafka producer.
       :host - IP address of kafka cluster."
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
        agent-message (dissoc agent-message :blobs)
        stream (java.io.ByteArrayOutputStream.)]
    (write-chunk! "JSON" (.getBytes (json/generate-string agent-message)) stream)
    (doseq [blob blobs]
      (write-chunk! "BLOB" blob stream))
    (.toByteArray stream)))

(defn message-with-optional-blobs
  "Associate blobs to agent message if they exist."
  [agent-message blobs]
  (if (empty? blobs)
    agent-message
    (assoc agent-message :blobs blobs)))

(defn read-json-chunk
  [^ChunkWriter chunk]
  (let [body (String. (.-body chunk) "UTF-8")]
    (json/parse-string body true)))

(defn decode-json-chunk
  [previous-message ^ChunkWriter chunk]
  (if (empty? previous-message) (read-json-chunk chunk) previous-message))

(defn read-chunks!
  "Read an Agent message map with an optional :blobs list from a stream of chunks."
  [^java.io.InputStream payload-stream]
  (loop [chunks (ChunkReader/readAll payload-stream false)
         agent-message {}
         blobs []]
    (let [^ChunkWriter chunk (first chunks)
          more (rest chunks)]
      (case (if chunk (.-chunkType chunk) :done)
        :done
        (message-with-optional-blobs agent-message blobs)
        "JSON"
        (recur more (decode-json-chunk agent-message chunk) blobs)
        "BLOB"
        (recur more agent-message (conj blobs (.-body chunk)))
        (recur more agent-message blobs)))))

(defn deserialize-from-chunks
  "Deserialize an Agent message map with an optional :blobs list from a Kafka
  message payload in chunk format."
  [^"[B" payload-bytes]
  (read-chunks! (java.io.ByteArrayInputStream. payload-bytes)))

(defn boot-producer
  "Instantiate a kafka producer with the given config."
  [config]
  (kafka/producer
   (producer-config config)
   (kafka/keyword-serializer)
   (byte-array-serializer)))

(defn send!
  "Send a message using a kafka producer on the given topic."
  [producer topic message]
  (kafka/send!
   producer
   {:topic topic
    :value (serialize-to-chunks message)}))

(defn consumer-config
  "Generate the configuration for a kafka consumer.
     config - a map containing options for the kafka consumer.
       :host - IP address of kafka cluster.
       :group-id - the kafka group id of this consumer."
  [config]
  {:bootstrap.servers (:host config)
   :enable.auto.commit "true"
   :auto.commit.interval.ms "1000"
   :group.id (get config :group-id "flow")
   :auto.offset.reset "latest"})

(defn boot-consumer
  "Instantiate a kafka consumer with the given configuration."
  [config]
  (kafka/consumer
   (consumer-config config)
   (kafka/keyword-deserializer)
   (byte-array-deserializer)))

(defn handle-message
  "Handle a raw message from kafka, invoke the user provided handler on it then emit it to the
   browser over websockets.
     state - the overall state of the system.
     bus - a manifold bus that emits messages over websockets to the browser or client.
     producer - a kafka producer in case the user provided handler wishes to send a message.
     handle - a function"
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
                  agent-message (dissoc agent-message :blobs)
                  topic-message-pair {topic agent-message}]
              (log/info (:event agent-message "") topic-message-pair blob-note)
              (handle state bus producer topic agent-message)
              (swap! state assoc :last-message topic-message-pair)
              (bus/publish! bus topic (json/generate-string topic-message-pair)))))))
    (catch Exception e
      (log/error e))))

(defn poll!
  "Wait for message from kafka.
     consumer - a reference to a kafka consumer."
  [consumer]
  (try
    (kafka/poll! consumer poll-interval)
    (catch Exception e
      (log/error e))))

(defn consume
  "Poll for messages and send them to the handler when received.
     consumer - a reference to a kafka consumer.
     handle - function to call when a message is received."
  [consumer handle]
  (binding [factory/*json-factory* non-numeric-factory]
    ; TODO(jerry): Revisit unpacking the poll! result.
    (loop [records (poll! consumer)]
      (when-not (empty? records)
        (doseq [record records]
          (handle record)))
      (recur (poll! consumer)))))

(defn boot-kafka
  "Boot a kafka consumer and producer and also a message bus for communicating with the browser."
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
  "Handle messages from the browser or websocket client.
     state - overall state of the system.
     conn - connection to send messages to the browser.
     message - message to send to the browser."
  [state conn message]
  (condp = (:event message)
    "INITIALIZE"
    (let [last-message (get @(:state state) :last-message)]
      (stream/put! conn (json/generate-string last-message)))

    (let [producer (:producer state)]
      (send! producer "flow" message))))

(defn boot
  "Boot the system using the provided config map.
     config - map containing values needed to instantiate kafka and websockets."
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
  "Render the index.html page on request."
  [state]
  (fn [request]
    {:status 200
     :headers {"content-type" "text/html"}
     :body (slurp "resources/public/index.html")}))

(defn connect-websocket
  "Establish a connection to the websocket client."
  [request]
  (defer/catch
      (http/websocket-connection request)
      (fn [_])))

(def non-websocket
  {:status 400
   :headers {"content-type" "application/text"}
   :body "must connect using websocket request"})

(defn websocket-handler
  "Connect the websocket client to the kafka streams."
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
  "A set of basic routes that any browser communication requires."
  [state]
  [["/" :index (#'index-handler state)]
   ["/ws" :websocket (#'websocket-handler state)]])

(defn start
  "Start the system with the provided state."
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
