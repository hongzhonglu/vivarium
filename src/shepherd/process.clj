(ns shepherd.process
  (:require
   [clojure.java.io :as io]
   [taoensso.timbre :as log]))

(defn launch!
  "Launches an external process with the given command line and configuration.
     command - a sequence of strings representing the command to be invoked.
     config - a map containing various options that could be provided to the system.
       :dir - the directory to invoke the command in.
       :env - any environment variables that should be present for the invoked command.
       :clear - a flag to trigger the system to clear the environment before invocation."
  [command {:keys [dir env clear] :as config}]
  (let [builder (ProcessBuilder. (into-array String command))
        environment (.environment builder)]
    (log/info "launch" command "dir" dir)
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
  "Determine whether the given process is currently running or not.
     process - a map containing a reference to the process and also its associated IO streams."
  [process]
  (.isAlive (:process process)))

(defn wait
  "Wait for the given process to exit, until the provided timeout.
     process - a map containing a reference to the process and also its associated IO streams.
     timeout - how long to wait for the process before returning (default 30 seconds)."
  ([process] (wait process 30))
  ([process timeout]
   (.waitFor (:process process) 30 java.util.concurrent.TimeUnit/SECONDS)))

(defn kill!
  "Kill the given process.
     process - a map containing a reference to the process and also its associated IO streams."
  [process]
  (.destroy (:process process)))

(defn ensure-termination!
  "Wait for the process to terminate, then kill it if it has not exited.
     process - a map containing a reference to the process and also its associated IO streams.
     timeout - how long to wait for the process before returning (default 30 seconds)."
  ([process] (ensure-termination! process 30))
  ([process timeout]
   (wait process timeout)
   (when (alive? process)
     (log/info "process: had to kill" process)
     (kill! process))))

(defn stream-to
  "Pipe an output stream to a given input stream.
     process - a map containing a reference to the process and also its associated IO streams.
     from - the key for which of this process' output streams to pipe.
     to - an input stream to pipe to."
  ([process from to] (stream-to process from to {}))
  ([process from to args]
   (apply io/copy (process from) to args)))

(defn stream-string
  "Stream a process' associated output stream to a string.
     process - a map containing a reference to the process and also its associated IO streams.
     from - which output stream to emit (:out or :err, default :out)"
  ([process] (stream-string process :out))
  ([process from]
   (with-open [writer (java.io.StringWriter.)]
     (stream-to process from writer)
     (str writer))))

(defn stream-out
  "Stream a process' associated output stream to stdout.
     process - a map containing a reference to the process and also its associated IO streams.
     from - which output stream to emit (:out or :err)"
  [process from]
  (future (stream-to process from (System/out))))
