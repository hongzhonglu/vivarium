(ns shepherd.process
  (:require
   [clojure.java.io :as io]))

(defn launch
  [command {:keys [dir env clear]}]
  (let [builder (ProcessBuilder. (into-array String command))
        environment (.environment builder)]
    (when dir
      (.directory builder (io/file dir)))
    (when clear
      (.clear environment))
    (doseq [[k v] env]
      (.put environment k v))
    (let [process (.start builder)]
      {:in (.getOutputStream process)
       :out (.getInputStream process)
       :err (.getErrorStream process)
       :process process})))

(defn kill
  [process]
  (.destroy (:process process)))

(defn stream-out
  [process]
  (future (io/copy (:out process) (System/out))))
