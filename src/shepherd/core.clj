(ns shepherd.core
  (:require
   [manifold.deferred :as defer]
   [taoensso.timbre :as log]
   [cheshire.core :as json]
   [shepherd.message :as message]
   [shepherd.process :as process]
   [shepherd.agent :as agent]))

(defn handle-message
  "Handle messages from kafka.
     state - state of the system.
     node - a channel to send messages to the websocket client.
     nexus - a reference to the kafka cluster.
     topic - a string naming the topic this message came from.
     message - the actual message."
  [state node nexus topic message]
  (condp = (:event message)
    "ADD_AGENT"
    (agent/add-agent! state node nexus message)
    "REMOVE_AGENT"
    (agent/remove-agent! state node nexus message)
    "TRIGGER_ALL"
    (agent/control-agents! state node nexus "TRIGGER_AGENT" message)
    "PAUSE_ALL"
    (agent/control-agents! state node nexus "PAUSE_AGENT" message)
    "SHUTDOWN_ALL"
    (agent/control-agents! state node nexus "SHUTDOWN_AGENT" message)))

(defn view-agent
  "Extract the values from an agent map that we want to transmit to the websocket client."
  [agent]
  (let [config (select-keys agent [:agent_id :agent_type :agent_config])
        alive? (process/alive? (:agent agent))
        view (assoc config :alive alive?)]
    view))

(defn status-handler
  "Render the status of all agents present in the system as json."
  [state]
  (fn [request]
    (let [agents @(:agents state)
          status (mapv view-agent (vals agents))]
      {:status 200
       :headers {"content-type" "application/json"}
       :body (json/generate-string status)})))

(defn shepherd-routes
  "Any routes the shepherd requires."
  [state]
  [["/status" :status (#'status-handler state)]])

(defn boot
  "Boot the shepherd with the given configuration."
  [config]
  (let [agents (atom {})
        config (assoc config :routes shepherd-routes)
        state {:agents agents :config config}
        handle (partial handle-message state)
        config (assoc-in config [:kafka :handle-message] handle)
        flow (message/boot config)]
    (merge state flow)))

(defn -main
  [& args]
  (let [config (message/read-config "config/config.clj")
        {:keys [port kafka]} (:shepherd config)
        config (assoc config :port port)
        config (update config :kafka merge kafka)
        state (boot config)]
    (log/info "shepherd starting" (:port config))
    (log/info state)
    (message/start state)))
