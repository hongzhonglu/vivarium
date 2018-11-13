(ns shepherd.core
  (:require
   [manifold.deferred :as defer]
   [taoensso.timbre :as log]
   [flow.core :as flow]
   [shepherd.agent :as agent]))

(defn handle-message
  [state topic message]
  (condp = (:event message)
    "ADD_AGENT"
    (agent/add-agent state message)
    "REMOVE_AGENT"
    (agent/remove-agent state message)))

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
