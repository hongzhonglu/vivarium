(ns shepherd.core
  (:require
   [clojure.edn :as edn]
   [clojure.java.shell :as sh]
   [clojure.java.io :as io]
   [taoensso.timbre :as log]
   [kinsky.client :as kafka]
   [cheshire.core :as json]
   [cheshire.factory :as factory]
   [manifold.deferred :as defer]
   [shepherd.process :as process]))

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

(defn boot-consumer
  [config]
  (kafka/consumer
   (consumer-config config)
   (kafka/keyword-deserializer)
   (kafka/json-deserializer)))

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

(defn handle-message
  [state record]
  (try
    (if (= (first record) :by-topic)
      (let [topics (last record)]
        (doseq [[topic messages] topics]
          (doseq [message messages]
            (log/info topic ":" message)
            (let [value {topic (:value message)}]
              (log/info value))))))
    (catch Exception e
      (log/error (.getMessage e)))))

(defn boot-kafka
  [state config]
  (let [producer (boot-producer config)
        consumer (boot-consumer config)]
    (doseq [topic (:topics config)]
      (log/info "subscribing to topic" topic)
      (kafka/subscribe! consumer topic))
    {:producer producer
     :consumer (defer/future (consume consumer (partial handle-message state)))}))

(defn boot
  [{:keys [kafka] :as config}]
  (let [state (atom {:agents {}})
        {:keys [producer consumer]} (boot-kafka state kafka)]
    {:config config
     :producer producer
     :consumer consumer
     :state state}))

(defn execute-agent
  [path command options])
