{:port 41114
 :launch
 {:dir "/home/spanglry/Code/wcEcoli"
  :boot "agent.boot"}
 :kafka
 {:host "localhost:9092"
  :group-id "shepherd"
  :event-topic "shepherd-receive"
  :topics ["shepherd-receive"]}}
