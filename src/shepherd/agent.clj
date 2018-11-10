(ns shepherd.agent
  (:require
   [cheshire.core :as json]
   [taoensso.timbre :as log]
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

