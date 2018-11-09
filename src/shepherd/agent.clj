(ns shepherd.agent
  (:require
   [cheshire.core :as json]
   [shepherd.process :as process]))

(defn launch-agent
  [spec dir]
  (process/launch
   ["python" "-m" "environment.boot"
    "--id" (:id spec)
    "--type" (:type spec)
    "--config" (json/generate-string (:config spec))]
   {:dir dir}))

