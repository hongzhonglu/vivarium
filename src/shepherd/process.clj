(ns shepherd.process
  (:require
   [clojure.java.io :as io]
   [taoensso.timbre :as log]))

(defn launch!
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

(defn alive?
  [process]
  (.isAlive (:process process)))

(defn wait
  ([process] (wait process 30))
  ([process timeout]
   (future
     (.waitFor (:process process) 30 java.util.concurrent.TimeUnit/SECONDS))))

(defn kill!
  [process]
  (.destroy (:process process)))

(defn ensure-termination!
  ([process] (ensure-termination! process 30))
  ([process timeout]
   (wait (:process process) timeout)
   (when (alive? (:process process))
     (log/info "process: had to kill" process)
     (kill! (:process process)))))

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
  (future (stream-to process from (System/out))))
