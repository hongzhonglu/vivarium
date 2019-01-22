(ns shepherd.lens
  (:require
   [manifold.deferred :as defer]
   [manifold.stream :as stream]
   [taoensso.timbre :as log]
   [cheshire.core :as json]
   [shepherd.message :as message]
   [shepherd.process :as process]
   [shepherd.agent :as agent]))

(defn handle-message
  [state node nexus topic message])

(defn handle-client
  [state conn message]
  (condp = (:event message)
    "VISUALIZATION_INITIALIZE"
    (let [last-message (get @(:state state) :last-message)]
      (stream/put! conn (json/generate-string last-message)))

    "DIVIDE_CELL"
    (let [producer (:producer state)]
      (message/send! producer "cell-receive" message))

    (let [producer (:producer state)]
      (message/send! producer "shepherd-receive" message))))

(defn boot
  [config]
  (let [state {:config config}
        handle (partial handle-message state)
        config (assoc-in config [:kafka :handle-message] handle)
        config (assoc config :handle-client handle-client)
        flow (message/boot config)]
    (merge state flow)))

(defn -main
  [& args]
  (let [config (message/read-config "config/config.clj")
        {:keys [port kafka]} (:lens config)
        config (assoc config :port port)
        config (update config :kafka merge kafka)
        state (boot config)]
    (log/info config)
    (log/info "lens started on port" (:port config))
    (message/start state)))
