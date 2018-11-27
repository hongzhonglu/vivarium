{:port 41114
 :launch
 {:dir "/home/spanglry/Code/wcEcoli"
  :boot "environment.boot"}
 :kafka
 {:host "localhost:9092"
  :group-id "shepherd"
  :send "shepherd-receive"
  :subscribe ["shepherd-receive"]
  :topics
  {:shepherd_receive "shepherd-receive"
   :agent_receive "agent-receive"
   :environment_receive "environment-receive"
   :cell_receive "cell-receive"
   :visualization_receive "environment-state"}}}
