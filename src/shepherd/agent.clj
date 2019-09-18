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

(def default-invocation
  ["python" "-u" "-m"])

(defn assemble-boot
  "If boot is a string, append it to the default python invocation.
   If it is a seq, use it as the invocation (this allows boot to be any arbitrary program)."
  [boot]
  (if (string? boot)
    (conj default-invocation boot)
    boot))

(defn assemble-args
  [id type config]
  ["--id" id
   "--type" type
   "--config" config])

(defn launch-agent!
  "Launch a python process based on the given agent spec.
     spec - a map containing any information necessary to boot the agent.
       :agent_id - a unique id for the new agent.
       :agent_type - the type of the new agent to be invoked.
       :agent_config - a configuration map containing any values needed by the python
          script to boot the new agent.
     config - system configuration.
       :boot - either a string containing a python module to run, or a list containing
          the initial components of the invocation."
  [spec config]
  (let [agent-config (:agent_config spec)
        boot-arg (or
                  (:boot agent-config)
                  (get config :boot "agent.boot"))
        boot (assemble-boot boot-arg)
        serial (json/generate-string agent-config)
        args (assemble-args (:agent_id spec) (:agent_type spec) serial)
        command (mapv identity (concat boot args))]
    (process/launch! command config)))

(defn ensure-kafka-config
  "Enforce the presence of a kafka config if not present based on this system's kafka config."
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
          bare (dissoc message :blobs)]
      (if (pos? (count files))
        (assoc-in bare [:agent_config :files] files)
        bare))
    message))

(defn add-agent!
  "Instantiate a new agent with the given configuration.
     state - state of the system.
     node - a channel to send messages to the websocket client.
     nexus - a reference to the kafka cluster.
     message - the message containing any information necessary to boot the new agent.
       :agent_id - a unique id for the new agent.
       :agent_type - the type of the new agent to be invoked.
       :agent_config - a configuration map containing any values needed by the python script
          to boot the new agent."
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
  "Shutdown an agent with the given id.
     state - state of the system.
     node - a channel to send messages to the websocket client.
     nexus - a reference to the kafka cluster.
     id - id of the agent to shutdown."
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
  "Determine if a string matches the given prefix.
     prefix - a string representing the prefix we are looking for.
     s - the string to test if the prefix is present."
  [prefix s]
  (= prefix (.substring s 0 (count prefix))))

(defn remove-prefix!
  "Remove any agents with the given prefix.
     state - state of the system.
     node - a channel to send messages to the websocket client.
     nexus - a reference to the kafka cluster.
     prefix - prefix to hunt for agents to remove."
  [state node nexus prefix]
  (let [agent-ids (keys @(:agents state))
        matching (filter (partial match-prefix prefix) agent-ids)]
    (doseq [match matching]
      (shutdown-agent! state node nexus match))))

(defn remove-agent!
  "Remove agents based on the given message. If :agent_id is present, remove only that agent.
   If :prefix is present, remove all agents with the given prefix.
     state - state of the system.
     node - a channel to send messages to the websocket client.
     nexus - a reference to the kafka cluster.
     message - information to direct the removal of agents.
       :prefix - if present, remove all agents with this prefix.
       :agent_id - if present, remove the agent with the given id."
  [state node nexus message]
  (if-let [prefix (:prefix message)]
    (remove-prefix! state node nexus prefix)
    (shutdown-agent! state node nexus (:agent_id message))))

(defn control-agents!
  "Send a message to all agents in this shepherd.
     state - state of the system.
     node - a channel to send messages to the websocket client.
     nexus - a reference to the kafka cluster.
     event - which event this message represents.
     message - the message to send to all agents."
  [state node nexus event message]
  (let [topic (get-in state [:config :kafka :topics :agent_receive])]
    (doseq [[id agent] @(:agents state)]
      (message/send!
       nexus topic
       {:event event
        :agent_id id}))))
