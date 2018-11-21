{:port 41114
 :launch
 {:dir "/Users/rspangler/Code/wcEcoli"
  :boot "environment.boot"}
 :kafka
 {:host "localhost:9092"
  :group-id "shepherd"
  :send "shepherd-receive"
  :subscribe ["shepherd-receive"]
  :topics
  {:shepherd-receive "shepherd-receive"
   :agent-receive "agent-receive"}}}
