(ns shepherd.core
  (:require
   [manifold.deferred :as defer]
   [taoensso.timbre :as log]
   [flow.core :as flow]
   [shepherd.agent :as agent]))

(defn handle-message
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

(defn boot
  [config]
  (let [agents (atom {})
        state {:agents agents :config config}
        handle (partial handle-message state)
        config (assoc-in config [:kafka :handle-message] handle)
        flow (flow/boot config)]
    (assoc state :flow flow)))

(defn -main
  [& args]
  (let [config (flow/read-config "config/config.clj")
        state (boot config)]
    (log/info "shepherd starting" (:port config))
    (log/info state)
    (flow/start state)))
