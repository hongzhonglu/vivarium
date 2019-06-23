# Lens

manage a colony of distributed whole-cell agents

## Usage

1. Install and start your Zookeeper and Kafka servers per instructions in [CovertLab/wcEcoli/agent/README.md](https://github.com/CovertLab/wcEcoli/tree/master/agent/README.md).

2. Start the Shepherd server to manage agent processes:

   `> lein run`

   * You can check its status by opening a web browser onto [http://localhost:41114/status](http://localhost:41114/status).

3. Start the Lens visualization server:

   `> lein run -m shepherd.lens`

4. Open a web browser onto Lens [http://localhost:33332](http://localhost:33332).

   * The Run/Pause button controls the simulation clock.

   * The pop-up menu picks which molecular concentration field to view. "GLC[p]" and "CARBON-MONOXIDE[p]" are the most interesting.

5. Create agents per the instructions in [CovertLab/wcEcoli/environment/README.md](https://github.com/CovertLab/wcEcoli/tree/master/environment/README.md).


## License

Copyright © 2018 Covert Lab

Distributed under the Eclipse Public License either version 1.0 or (at
your option) any later version.
