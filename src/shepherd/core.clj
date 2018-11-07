(ns shepherd.core
  (:require
   [kinsky.client :as kafka]
   [clojure.java.sh :as sh]))

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

(defn consumer-config
  [config]
  {:bootstrap.servers (:host config)
   :enable.auto.commit "true"
   :auto.commit.interval.ms "1000"
   :group.id (get config :group-id "shepherd")
   :auto.offset.reset "latest"})

(defn producer-config
  [config]
  {:bootstrap.servers (:host config)})

(defn boot-producer
  [config]
  (kafka/producer
   (producer-config config)
   (kafka/keyword-serializer)
   (kafka/json-serializer)))

(defn send!
  [producer topic message]
  (kafka/send!
   producer
   {:topic topic
    :value message}))

(defn consume
  [consumer handle]
  (binding [factory/*json-factory* non-numeric-factory]
    (loop [records
           (try
             (kafka/poll! consumer poll-interval)
             (catch Exception e
               (log/error (.getMessage e))))]
      (when-not (empty? records)
        (doseq [record records]
          (handle record)))
      (recur (kafka/poll! consumer poll-interval)))))

(defn boot-kafka
  [state config]
  (let [bus (bus/event-bus)
        producer (boot-producer config)
        consumer (boot-consumer config)]
    (doseq [topic (:topics config)]
      (log/info "subscribing to topic" topic)
      (kafka/subscribe! consumer topic))
    {:bus bus
     :producer producer
     :consumer (defer/future (consume consumer (partial handle-message state bus)))}))

(defn boot
  [{:keys [kafka] :as config}]
  (let [state (atom {:last-message {}})
        {:keys [bus producer consumer]} (boot-kafka state kafka)]
    {:config config
     :bus bus
     :producer producer
     :consumer consumer
     :state state}))

(defn execute-agent
  [path ])
