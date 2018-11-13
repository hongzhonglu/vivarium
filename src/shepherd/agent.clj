(ns shepherd.agent
  (:require
   [cheshire.core :as json]
   [taoensso.timbre :as log]
   [flow.core :as flow]
   [shepherd.process :as process]))

(defn launch-agent
  [spec config]
  (let [serial (json/generate-string (:config spec))]
    (log/info serial)
    (process/launch
     ["python" "-u" "-m" (get config :boot "agent.boot")
      "--id" (:id spec)
      "--type" (:type spec)
      "--config" serial]
     config)))

(defn add-agent
  [state message]
  (log/info "add agent:" message state)
  (let [record (select-keys message [:id :type :config])
        launch-config (get-in state [:config :launch])
        _ (log/info "launch config" launch-config)
        born (launch-agent message launch-config)
        record (assoc record :agent born)]
    (swap! (:agents state) assoc (:id record) record)
    (process/stream-out (:agent record) :out)
    (process/stream-out (:agent record) :err)))

(declare remove-prefix)

(defn remove-agent
  [state message]
  (if (:prefix message)
    (remove-prefix state message)
    (let [id (:id message)
          agent (get @(:agents state) id)
          producer (get-in state [:flow :producer])
          topic (get-in state [:config :kafka :topics :agent-receive])]
      (flow/send! producer topic {:event "SHUTDOWN_AGENT" :agent_id id}))))
