(ns shepherd.process
  (:require
   [clojure.java.io :as io]
   [taoensso.timbre :as log]))

(defn launch
  [command {:keys [dir env clear] :as config}]
  (let [builder (ProcessBuilder. (into-array String command))
        environment (.environment builder)]
    (.redirectErrorStream builder true)
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

(defn stream-to
  ([process from to] (stream-to process from to {}))
  ([process from to args]
   (apply io/copy (process from) to args)))

(defn stream-string
  ([process] (stream-string process :out))
  ([process source]
   (with-open [writer (java.io.StringWriter.)]
     (stream-to process source writer)
     (str writer))))

(defn stream-out
  [process from]
  (stream-to process from (System/out)))
