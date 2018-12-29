(defproject shepherd "0.0.1"
  :description "manage a flock of distributed agents"
  :url "https://github.com/CovertLab/shepherd"
  :license {:name "Eclipse Public License"
            :url "http://www.eclipse.org/legal/epl-v10.html"}
  :dependencies [[org.clojure/clojure "1.8.0"]
                 [org.slf4j/slf4j-simple "1.7.25"]
                 [com.taoensso/timbre "4.8.0"]
                 [aleph "0.4.6"]
                 [clj-http "3.7.0"]
                 [polaris "0.0.19"]
                 [spootnik/kinsky "0.1.22"]]
  :main shepherd.core
  :java-source-paths ["java/src"])
