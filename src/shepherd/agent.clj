(ns shepherd.agent
  (:require
   [cheshire.core :as json]
   [taoensso.timbre :as log]
   [flow.core :as flow]
   [shepherd.process :as process]))

(defn launch-agent!
  [spec config]
  (let [serial (json/generate-string (:config spec))]
    (log/info serial)
    (process/launch!
     ["python" "-u" "-m" (get config :boot "agent.boot")
      "--id" (:id spec)
      "--type" (:type spec)
      "--config" serial]
     config)))

(defn add-agent!
  [state message]
  (log/info "add agent:" message state)
  (let [record (select-keys message [:id :type :config])
        launch-config (get-in state [:config :launch])
        _ (log/info "launch config" launch-config)
        born (launch-agent! message launch-config)
        record (assoc record :agent born)]
    (swap! (:agents state) assoc (:id record) record)
    (process/stream-out (:agent record) :out)
    (process/stream-out (:agent record) :err)))

(defn shutdown-agent!
  [state id]
  (let [agent (get @(:agents state) id)
        producer (get-in state [:flow :producer])
        topic (get-in state [:config :kafka :topics :agent-receive])]
    (flow/send!
     producer topic
     {:event "SHUTDOWN_AGENT"
      :agent_id id})
    (process/ensure-termination! (:agent agent))
    (swap! (:agents state) dissoc id)))

(defn match-prefix
  [prefix s]
  (= prefix (.substring s 0 (count prefix))))

(defn remove-prefix!
  [state prefix]
  (let [agent-ids (keys @(:agents state))
        matching (filter (partial match-prefix prefix) agent-ids)]
    (doseq [match matching]
      (shutdown-agent! state match))))

(defn remove-agent!
  [state message]
  (if-let [prefix (:prefix message)]
    (remove-prefix! state prefix)
    (shutdown-agent! state (:id message))))

(defn control-agents!
  [state message]
  (let [producer (get-in state [:flow :producer])
        topic (get-in state [:config :kafka :topics :agent-receive])]
    (doseq [agent @(:agents state)]
      (flow/send!
       producer topic
       {:event (:event message)
        :agent_id (:id agent)}))))
