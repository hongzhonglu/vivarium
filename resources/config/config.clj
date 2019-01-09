{:lens
 {:port 33332
  :kafka
  {:group-id "lens-visualization"
   :event-topic "environment-state"
   :subscribe ["environment-state"]}}

 :shepherd
 {:port 41114
  :kafka
  {:group-id "shepherd"
   :event-topic "shepherd-receive"
   :subscribe ["shepherd-receive"]}}

 :launch
 {:dir "../wcEcoli"
  :boot "environment.boot"}

 :kafka
 {:host "localhost:9092"
  :topics
  {:shepherd_receive "shepherd-receive"
   :agent_receive "agent-receive"
   :environment_receive "environment-receive"
   :cell_receive "cell-receive"
   :visualization_receive "environment-state"}}}
