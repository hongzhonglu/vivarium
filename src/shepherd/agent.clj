(ns shepherd.agent
  (:require
   [cheshire.core :as json]
   [taoensso.timbre :as log]
   [shepherd.message :as message]
   [shepherd.process :as process])
  (:import
    [java.io IOException]
    [java.nio.file Files OpenOption]
    [java.nio.file.attribute FileAttribute]))

(defn launch-agent!
  [spec config]
  (let [serial (json/generate-string (:agent_config spec))]
    (process/launch!
     ["python" "-u"
      "-m" (get config :boot "agent.boot")
      "--id" (:agent_id spec)
      "--type" (:agent_type spec)
      "--config" serial]
     config)))

(defn ensure-kafka-config
  [state message]
  (if (get-in message [:agent_config :kafka_config])
    message
    (let [all-kafka (get-in state [:config :kafka])
          kafka (select-keys all-kafka [:host :topics])
          kafka (assoc kafka :subscribe [])]
      (assoc-in message [:agent_config :kafka_config] kafka))))

(defn write-temp-blob!
  "Write a byte array to a new temp file and return its file path."
  [^"[B" blob]
  (try
    (let [path (Files/createTempFile "agent" ".blob" (make-array FileAttribute 0))]
      (.deleteOnExit (.toFile path))
      (when blob
        (Files/write path blob (make-array OpenOption 0)))
      (.toString path))
    (catch IOException e (log/error e))
    (catch SecurityException e (log/error e))))

(defn blobs-to-temp-files!
  "Move message's :blobs to temp :files as positional args."
  [message]
  (if-let [blobs (:blobs message)]
    (let [files (mapv write-temp-blob! blobs)
          msg (dissoc message :blobs)]
      (if (pos? (count files))
        (assoc-in msg [:agent_config :files] files)
        msg))
    message))

(defn add-agent!
  [state node nexus message]
  (let [record (select-keys message [:agent_id :agent_type :agent_config])
        message (blobs-to-temp-files! (ensure-kafka-config state message))
        launch-config (get-in state [:config :launch])
        born (launch-agent! message launch-config)
        record (assoc record :agent born)]
    (swap! (:agents state) assoc (:agent_id record) record)
    (process/stream-out (:agent record) :out)
    (process/stream-out (:agent record) :err)))

(defn shutdown-agent!
  [state node nexus id]
  (let [agent (get @(:agents state) id)
        topic (get-in state [:config :kafka :topics :agent_receive])]
    (message/send!
     nexus topic
     {:event "SHUTDOWN_AGENT"
      :agent_id id})
    (process/ensure-termination! (:agent agent))
    (swap! (:agents state) dissoc id)))

(defn match-prefix
  [prefix s]
  (= prefix (.substring s 0 (count prefix))))

(defn remove-prefix!
  [state node nexus prefix]
  (let [agent-ids (keys @(:agents state))
        matching (filter (partial match-prefix prefix) agent-ids)]
    (doseq [match matching]
      (shutdown-agent! state node nexus match))))

(defn remove-agent!
  [state node nexus message]
  (if-let [prefix (:prefix message)]
    (remove-prefix! state node nexus prefix)
    (shutdown-agent! state node nexus (:agent_id message))))

(defn control-agents!
  [state node nexus event message]
  (let [topic (get-in state [:config :kafka :topics :agent_receive])]
    (doseq [[id agent] @(:agents state)]
      (message/send!
       nexus topic
       {:event event
        :agent_id id}))))
