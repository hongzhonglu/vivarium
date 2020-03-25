===============
Getting Started
===============

-------------------------
Download and Installation
-------------------------

This guide assumes that you have access to an administrator or sudoer
account on a macOS or Linux system.

Getting Organized
=================

Creating Enclosing Directory
----------------------------

Create a ``vivarium_work`` folder anywhere you like. But for installing
some third-party software, everything we do will occur inside this
folder.

.. _pythonpath:

Setting PYTHONPATH
------------------

Vivarium needs the root of the repository to be in your ``PYTHONPATH``
environment variable so that Python can find Vivarium. To make this easy
to set, we suggest adding this line to your shell startup file:

.. code-block:: bash

    alias pycd='export PYTHONPATH="$PWD:$PYTHONPATH"'

Now when you are about to work on Vivarium, navigate to the root of the
Vivarium repository (``vivarium_work/vivarium``) and run ``pycd`` in
your terminal. You will need to do this for each terminal window you
use.

Installing Dependencies
=======================

Below, we list the dependencies Vivarium requires, how to check whether
you have them, how to install them, and in some cases, how to set them
up for Vivarium.  Make sure you have each of them installed.

Python 3
--------

The Vivarium code is written in Python.

*Check Installation*

.. code-block:: console

    $ python --version
    Python <version>

Make sure you see a version at least 3.6.

*Install*

Download the latest installer from https://www.python.org/downloads/

**Open JDK 8**

Zookeeper and Kafka, which we will address later, require that you have
Java installed.

*Check Installation*

.. code-block:: console

    $ java -version
    java version <version> <date>

Ensure the version is at least 1.8.

*Install*

Download the latest JDK installer from
https://www.oracle.com/java/technologies/javase-downloads.html.  As of
writing the latest is Java SE 14. Download and run the appropriate
installer for your platform. Then you need to set the ``JAVA_HOME``
environment variable, for instance by adding ``export
JAVA_HOME=$(/usr/libexec/java_home)`` to your startup shell file (e.g.
``~/.bash_profile`` or ``~/.profile``).

Zookeeper and Kafka
-------------------

Kafka is a message passing system that allows decoupling of message
senders and message receivers. It does this by providing two
abstractions, a Consumer and a Producer. A Consumer can subscribe to any
number of "topics" it will receive messages on, and a Producer can send
to any topics it wishes. Topics are communication "channels" between
processes that otherwise do not need to know who is sending and
receiving these messages. Vivarium uses Kafka to pass messages between
actors, for example between a cell and its environment.

Kafka is built on Zookeeper, which is a service for synchronizing access
to a hierarchy of key/value pairs called "nodes". The Kafka and
Zookeeper services are meant to be deployed in multiple instances in a
server cluster, but there is a binary distribution that you can run
locally for development.

*Check Installation*

For this guide, we will not install Kafka globally on your system.
Instead, we will store the Kafka program in ``vivarium_work`` and run
the executable directly. This means you almost certainly need to install
it, even if you use Kafka already.

*Install*

#. Download Kafka from https://kafka.apache.org/downloads, choosing the
   latest version. This will give you a ``.tgz`` archive file that
   includes both Kafka and Zookeeper.
#. Unarchive this file into ``vivarium_work`` to create a folder like
   ``vivarium_work/kafka_2.11-2.0.0/``. Your folder name will likely
   change slightly to match your version of Kafka.
#. Create a shell script ``vivarium_work/zookeeper.sh`` with the
   following content:

   .. code-block:: bash

        #!/bin/bash
        
        ./kafka_2.11-2.0.0/bin/zookeeper-server-start.sh \
            ./kafka_2.11-2.0.0/config/zookeeper.properties

#. Create a shell script ``vivarium_work/kafka.sh`` with the following
   content:

   .. code-block:: bash

        #!/bin/bash

        ./kafka_2.11-2.0.0/bin/kafka-server-start.sh \
            ./kafka_2.11-2.0.0/config/server.properties \
            --override listeners=PLAINTEXT://127.0.0.1:9092

   Overriding the "listeners" address like this allows connections to
   the Kafka server to withstand network DHCP address changes and the
   like.
#. Make the scripts executable like this:

   .. code-block:: console

        $ chmod 700 vivarium_work/kafka.sh
        $ chmod 700 vivarium_work/zookeeper.sh

   Now you can easily start and stop the zookeeper and Kafka servers
   like this:

   .. code-block:: console

        $ vivarium_work/zookeeper.sh
        $ vivarium_work/kafka.sh

   Note that Zookeeper must be started before Kafka. Also note that
   these two commands must be run in separate terminals. To shut them
   down, you can just use CTRL-C to kill the processes.

   .. note:: Make sure you shut down Kafka before Zookeeper!  If you
       shut down Zookeeper first, Kafka will refuse to quit. You can
       then force it to stop with ``kill -9``.
 
MongoDB
-------

We use a MongoDB database to store the data collected from running
simulations. This can be a a remote server, but for this guide we will
run a MongoDB server locally.

*Check Installation*

.. code-block:: console

    $ mongod --version
    db version v4.2.3
    ...

Make sure you see a version at least 4.2.

*Install*

If you are on macOS, you can install MongoDB using `Homebrew
<https://brew.sh>`_. You will need to add the MongoDB tap following the
instructions `here <https://github.com/mongodb/homebrew-brew>`_.

If you are on Linux, see the MongoDB documentation's `instructions
<https://docs.mongodb.com/manual/administration/install-on-linux/>`_.

*Setup*

There are any number of ways to get a MongoDB server up and running
locally. Here is one way:

#. Create a folder ``vivarium_work/mongodb``. This is where the database
   will be stored. We store the database here instead of at the default
   location in ``/usr/local/var/mongodb`` to avoid permissions issues if
   you are not running as an administrator.
#. Make a copy of the mongod configuration file so we can make changes:

   .. code-block:: console

      $ cp /usr/local/etc/mongod.conf vivarium_work/mongod.conf

   Note that your configuration file may be somewhere slightly
   different. Check the MongoDB documentation for your system.
#. In ``vivarium_work/mongod.conf`` change the path after ``dbPath:`` to
   point to ``vivarium_work/mongodb``.
#. Create a shell script ``vivarium_work/mongo.sh`` with the following
   content:

   .. code-block:: bash

      #!/bin/bash

      mongod --config mongodb.conf

#. Make the script executable:

   .. code-block:: console

        $ chmod 700 vivarium_work/mongo.sh

   Now you can launch MongoDB by running this script:

   .. code-block:: console

        $ vivarium_work/mongo.sh

.. todo:: Use python -m agent.boot --host ip.to.remote.cluster:9092
    for remote Kafka services

GNU Linear Programming Kit (GLPK)
---------------------------------

.. todo:: What is GLPK used for?

One of the Python packages we will install later, ``swiglpk``, requires
that GLPK already be installed on your system.

*Check Installation*

Currently we don't have a way to check whether ``glpk`` is installed. If
you think you already have it, you can proceed with the installation and
watch for an error about missing ``glpk``.

.. todo:: Check GLPK installation

*Install*

If you use Homebrew, you
can install GLPK like this:

.. code-block:: console

    $ brew install glpk

Otherwise, follow the installation instructions on the GLPK
`homepage <https://www.gnu.org/software/glpk>`_.

Leiningen
---------

Our simulation runs each cell on its own thread, and we use Leiningen
to manage these threads.

*Check Installation*

To check whether Leiningen is installed, run:

.. code-block:: console

    $ lein --version
    Leiningen <version> ...

You may also see a deprecation warning from Java HotSpot, which you can
ignore. Make sure the version is at least 2.9.

*Install*

To install Leiningen, follow the instructions on its `website
<https://leiningen.org/>`_. Alternatively, you can also install the
``leiningen`` formula on Homebrew.

Download and Setup Vivarium
===========================

Download the Code
-----------------

The Vivarium code is available on `GitHub
<https://github.com/CovertLab/vivarium>`_. Move into your
``vivarium_work`` directory and clone the repository to
download the code

.. code-block:: console
    
    $ cd vivarium_work
    $ git clone https://github.com/CovertLab/vivarium.git

This will create a ``vivarium`` folder inside ``vivarium_work``. All the
code for Vivarium is inside this ``vivarium`` folder.

Installing Python Packages
--------------------------

Above we installed all the non-Python dependencies, but we still have to
install the Python packages Vivarim uses.

#. Move into the ``vivarium`` folder created when you cloned the
   repository.
#. (optional) Create and activate a virtual environment:

   .. code-block:: console

        $ python3 -m venv venv
        ...
        $ source venv/bin/activate
#. Install packages

   .. code-block:: console

        $ pip install -r reqirements.txt

   If you encounter problems installing numpy and/or scipy, try this
   instead:

   .. code-block:: console

        $ pip install -r requirements.txt --no-binary numpy,scipy
        $ pip install numpy
        $ pip install scipy

Now you are all set to run Vivarium!

---------------
Run Simulations
---------------

Some Terminology: Processes and Composites
==========================================

In Vivarium, we break our cell models into *processes*. Each process
models part of the cell's function. For example, we have processes for
metabolism, transcription, and translation in Vivarium. We can combine
these processes into *composites* that model a cell with all of the
functionality modeled by the included processes. For example, we could
compose transcription and translation to create a more complete gene
expression model.

.. todo:: Link to topical guide on processes and composites

In Vivarium, we store processes in ``vivarium/vivarium/processes`` and
composites in ``vivarium/vivarium/composites``.

Running Processes and Composites in Isolation
=============================================

Every process and composite can be run by itself. While this is not
particularly useful for modeling cells or colonies, it may help
illustrate what the process or composite does. To run a process or
composite, you can simply execute the Python file that defines it. For
example, we can run the degradation process like this:

.. code-block:: console

    $ python vivarium_work/vivarium/vivarium/processes/degradation.py
    ...

.. note:: If you get errors from Python about being unable to find
    ``vivarium``, make sure you've set your PYTHONPATH correctly. See
    :ref:`pythonpath` for details.

Don't worry about the output--it's only useful for developers. You will
see that a new folder has been created at
``vivarium_work/vivarium/vivarium/out/tests``. This is where we store
the output from running processes and composites in isolation. For the
degradation process, the output is in the ``degradation`` folder inside
``tests``. Here you'll find a ``simulation.png`` file that looks like
this:

.. image:: ./_static/degradation_plots.png
    :width: 100%
    :alt: Four columns of plots, each of which has the plotted value on
        the y-axis and time on the x-axis. In the first column, we see
        the concentration of transcripts decreasing linearly with time,
        while in the third column concentrations of the four RNA
        nucleotides increase linearly with time. In the second column a
        plot of the concentration of endoRNAse is a horizontal line, and
        in the fourth column plots of metrics like density, volume, and
        mass are all constant.

This shows quickly that this degradation process removes transcripts and
returns the RNA nucleotides to the cell.

Some processes also produce the data from which the plots were produced.
This is saved in ``simulation_data.csv``. Try running the
``convenience_kinetics`` process to see how this works!

Lastly, try running the ``flagella_expression`` composite like this:

.. code-block:: console

    $ python vivarium_work/vivarium/vivarium/composites/flagella_expression.py

Now in the ``flagella_expression_composite`` in ``tests``, you should
see an image containing a plot like this:

.. image:: ./_static/flagella_expression_aa_plot.png
    :width: 100%
    :alt: Five plots showing the concentrations of various polymerases,
        nucleotides, amino acids, transcripts, and proteins over time.
        The amino acid plot shows one amino acid running out first.

Notice that even from this minimal simulation, we can tell which amino
acid is limiting! In this case the colors are so similar that it's hard
to tell, but the limiting amino acid is either alanine or leucine.

.. todo:: Is alanine or leucine limiting?

.. _agents-in-terminal-windows:

Running Agents in Terminal Windows
==================================

.. note:: Running agents separately in terminal windows is very useful
    for debugging because it lets you see the output from each agent.

Terminology: Agents
-------------------

Vivarium is heavily influenced by agent-based modeling, in which the
model consists of individual agents interacting with each other. In
Vivarium, each cell is an agent. The environment is also an agent. These
agents interact with each other by passing messages through Kafka.

.. todo:: Link to more comprehensive topical guide

How to Run Agents
-----------------

Each agent runs on its own thread. We do this because each agent can be
as complex as an entire whole-cell model, so the entire simulation
cannot be run on a single thread. Shepherd can manage these threads for
you; importantly, you must use Shepherd if your simulation will require
creating or deleting threads. Cell division, for example, involves
stopping the mother cell's thread and starting two new threads, one for
each daughter cell, so division requires Shepherd.

.. todo:: Link to using Shepherd

That said, you *can* run agents on your own instead of using Shepherd.

.. note:: If you run a simulation using this method that includes
    stopping and/or starting agents, the agents will stop, but new ones
    will not start. For example if your cell divides, the agent you
    started for the mother cell will stop, but the daughter cells will
    not be created.

We will run each agent in its own terminal window to mimic the threads
that Shepherd would create. Let's see how!

First we need to get all our servers running. Do each of the following
in a separate terminal window:

#. Start Zookeeper:

   .. code-block:: console
   
        $ vivarium_work/zookeeper.sh
        ...
        ... INFO binding to port 0.0.0.0/0.0.0.0:2181 ...

#. Start Kafka:

   .. code-block:: console
   
        $ vivarium_work/kafka.sh
        ...
        ... INFO [KafkaServer id=0] started (kafka.server.KafkaServer)

   You should also see som text print out on the Zookeeper window. You
   might see some ``NoNode`` warnings--these are safe to ignore.

   .. note:: Zookeeper must be started before Kafka!

#. Start MongoDB:

   .. code-block:: console
   
        $ vivarium_work/mongo.sh

   There shouldn't be any output.

   .. note:: Alternatively, if you installed MongoDB using Homebrew, you
       can tell Homebrew to always run a MongoDB server by running:

       .. code-block:: console

            $ brew services start mongodb/brew/mongodb-community

       Now a MongoDB server will be started automatically once you
       login. Then you can skip the step of starting MongoDB in the
       future.

Now we can create our agents. We create an agent like this:

.. code-block:: console

    $ python -m vivarium.environment.boot --type <type> --id <id> [--outer-id <outId>]

.. note:: If you get errors from Python about being unable to find
    ``vivarium``, make sure you've set your PYTHONPATH correctly. See
    :ref:`pythonpath` for details.

where ``<type>`` is the agent type, ``<id>`` is the identifier for this
agent, and ``<outId>`` is an optional argument that stipulates that the
agent should be placed inside the agent with identifier ``<outId>``.
This outer agent will almost always be an environment. There is also an
optional ``--config '{...}'`` argument you can use to configure the
agent.

.. todo:: Link to information on configuration

Ther are many possible agent types. To see them all check out the help
text like this:

.. code-block:: console

    $ python -m vivarium.environment.boot --help

.. todo:: Point to the autogenerated docs for the agents

Here's an example of running a simulation of a simple environment with
three cells that consume glucose and lactose. We will initialize the
environment with glucose and lactose, and as the glucose is depleted, we
should see the cells shift to consuming lactose.

.. WARNING:: This example doesn't work yet. The general process is
    correct, but the particular agent types are not.

.. todo:: Instructions for debugging in this mode

#. First, let's create a ``ecoli_core_glc`` environment agent. This is a
   kind of lattice environment. Lattice environments discretize the
   simulation space into a two-dimensional grid, each region of which
   has the same depth. Each region has uniform metabolite
   concentrations, but metabolite concentrations differ between regions,
   letting us model a continuous distribution of concentrations. A
   diffusion process in the environment tends to make the space
   homogeneous. We start this agent like this:

   .. code-block:: console
   
        $ python -m vivarium.environment.boot --type ecoli_core_glc --id env
        environment started

   .. note:: Wait for the ``environment started`` to show up before
       proceeding. Otherwise there won't be an environment to add the
       cells to!

#. Next, let's create three cell agents. These agents will be of type
   ``shifter`` because they will initially consume glucose, but when
   glucose concentrations drop, they will start consuming lactose. We
   create these agents like this:

   .. code-block:: console
   
      $ python -m vivarium.environment.boot --type shifter --id c1 --outer-id env
      $ python -m vivarium.environment.boot --type shifter --id c2 --outer-id env
      $ python -m vivarium.environment.boot --type shifter --id c3 --outer-id env

   After creating each cell agent, you should see in both the cell and
   the environment's terminal windows a message from the cell to the
   environment declaring itself:

   .. code-block:: console

        <-- environment-receive CELL_DECLARE [shifter c1]: {'event':
        'CELL_DECLARE', 'agent_id': 'env', 'inner_id': 'c1',
        'agent_config': { ... }, 'state': {'volume': 1.0}}

   And a message from the environment back to the cell:

   .. code-block:: console
   
        <-- cell-receive ENVIRONMENT_SYNCHRONIZE [glc_lct env]:
        {'event': 'ENVIRONMENT_SYNCHRONIZE', 'inner_id': 'c1',
        'outer_id': 'env', 'state': { ... }}

#. Now we can start the simulation!

   .. code-block:: console
   
        $ python -m vivarium.environment.control run --id env

   The simulation will stop on its own once the environment agent hits
   the end of its programmed timeline. However, you can pause, run,
   and shutdown the simulation like this as well:

   .. code-block:: console
   
        $ python -m vivarium.environment.control pause --id env
        $ python -m vivarium.environment.control run --id env
        $ python -m vivarium.environment.control shutdown

.. todo:: Fix this tutorial, as currently it fails.

.. todo:: Add example output from this tutorial.

Using Shepherd
==============

The usual way to start the simulation is to use Shepherd, which spawns
agents in new threads as requested via Kafka messages so you don't have
to launch each agent in its own terminal tab. Furthermore, this enables
cell division wherein a cell agent process ends and two new ones begin.
But to debug an agent, see the :ref:`agents-in-terminal-windows`
instructions above.

Let's take a look at an example of using Shepherd. We'll be able to
model cells dividing!

.. todo:: Reference composites in this and the previous tutorial

.. WARNING:: This tutorial has not yet been tested!

#. Launch Shepherd in a separate terminal window:

   .. code-block:: console

        $ lein run

#. For our environment, let's make a ``lattice`` agent:

   .. code-block:: console
   
        $ python -m vivarium.environment.boot --type lattice --id env
        environment started

   .. note:: Wait for the ``environment started`` to show up before
       proceeding. Otherwise there won't be an environment to add the
       cells to!

#. Next, let's create a cell agent of type ``growth_division``, which
   can grow and divide.

   .. code-block:: console
   
      $ python -m vivarium.environment.boot --type growth_division --id c --outer-id env

#. Now we can start the simulation!

   .. code-block:: console
   
        $ python -m vivarium.environment.control run --id env
