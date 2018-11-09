(ns shepherd.agent
  (:require
   [cheshire.core :as json]
   [shepherd.process :as process]))

(defn launch-agent
  [agent-spec]
  (process/launch
   ["python" "-m" "environment.boot"
    "--id" (:agent-id agent-spec)
    "--type" (:agent-type agent-spec)
    "--config" (json/generate-string (:agent-config agent-spec))]))

(defn add-agent
  [shepherd kafka agent-spec])
