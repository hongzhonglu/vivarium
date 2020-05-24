==========
Composites
==========

To model a cell faithfully, we need to simulate processes running
concurrently and interacting with each other. For example, we might want
to model the synthesis of ATP by a metabolism process concurrently with
the consumption of ATP by a transport process such as a sodium pump. In
Vivarium, we model this by creating a :term:`composite`.

--------------------
Processes and Stores
--------------------

A composite models cellular functions running concurrently. Each
function is modeled as a :term:`process` in the composite. To let these
processes interact, for example by producing and consuming a shared
resource like ATP, the processes share parts of the state, called
:term:`stores`. Each store is a collection of :term:`variables` such as
cytoplasmic ATP concentration, and in the composite we define which
process operates on which stores.

In our ATP example, we might assign a "cytoplasm" store to both the
metabolism and sodium pump processes. Now when the composite is
simulated, the metabolism and sodium pump processes will be changing the
same variable, the ATP concentration in the cytoplasm store. This means
that if metabolism is running more slowly, the cytoplasmic ATP
concentration variable will be lower, so the sodium pump will export
less sodium. Thus, shared stores in composites let us simulate
interacting concurrent processes.

.. note:: Each variable must belong to exactly one store. This means
    that stores are partitions of the total model state.

How Processes and Stores are Implemented
========================================

Processes
---------

Processes are written as classes that inherit from
:py:class:`vivarium.core.process.Process`.  To create a
composite, we create instances of these classes to make the processes we
want to compose. If a process is configurable, we might provide a
dictionary of configuration options. For information on configuring a
process, see the process's documentation.

.. todo:: Point to an example documented process and its usage

Each process in a composite is given a name, which is used in the
topology we describe below.

Stores
------

Each store is implemented as a Python dictionary where the keys are the
variable names and the values are the variable values. We also have
wrapper :py:class:`vivarium.core.process.Store` objects that hold
these dictionaries.

----------------------------
Ports Make Processes Modular
----------------------------

We don't want process creators to worry about what kind of composite
their process will be used in. Conversely, if you are creating a
composite, you should be able to use any processes you like, even if
they weren't written with your use case in mind. Vivarium achieves this
modularity with :term:`ports`.

Each process has a list of named ports, one for each store it is
expecting. The process can perform all of its computations in terms of
these ports, and the process also provides its update using port names.
This means that each process can be used with any collection of stores,
making processes modular.

This modularity is analogous to the modularity of Python functions. Each
process can be thought of as a function like this:

.. code-block:: python

    def sodium_pump(cytoplasm, extracellularSpace):
        ...
        return "Update: Decrease ATP concentration in cytoplasm by x mM"

A function's modularity comes from the fact that we can pass in different
objects for the ``cytoplasm`` parameter, even objects the function
authors hadn't thought of. ``cytoplasm`` is like the port, to which we
can provide any store we like.

How Processes Define Ports
==========================

A process specifies its port names in its constructor by calling the
superclass (:py:class:`vivarium.core.process.Process`)
constructor. For example, the
:py:class:`vivarium.processes.convenience_kinetics.ConvenienceKinetics`
class contains this line:

.. code-block:: python

    super(ConvenienceKinetics, self).__init__(ports, parameters)

The ``ports`` variable takes the form of a dictionary with port names as
keys and lists of variable names as values. For example, if ``ports``
looked like this:

.. code-block:: python

    {
        'cytoplasm': ['ATP', 'sodium'],
        'extracellular': ['sodium']
    }

then the process would be declaring that it cares about the ``ATP`` and
``sodium`` variables in the ``cytoplasm`` port and the ``sodium``
variable in the ``extracellular`` port. When the process is asked to
provide an update to the model state, it is only provided the variables
it specifies. For example, it might get a model state like this:

.. code-block:: python

    {
        'cytoplasm': {
            'ATP': 5.0,
            'sodium': 1e-2,
        },
        'extracellular': {
            'sodium': 1e-1,
        },
    }

This would happen even if the store linked to the ``cytoplasm`` port
contained more variables. We call this stripping-out of variables the
process doesn't need :term:`masking`.

----------
Topologies
----------

How do we specify which store goes with which port? To continue the
function analogy from above, we need something analogous to this:

.. code-block:: python

    cell = Cell()
    bloodVessel = BloodVessel()
    # We need something like the line below
    update = sodium_pump(cytoplasm=cell, extracellularSpace=bloodVessel)

When we call ``sodium_pump``, we specify which objects go with which
parameters. Analogously, we specify the mapping between ports and stores
using a :term:`topology`.

Defining Topologies
===================

Topologies are defined as dictionaries with process names as keys and
dictionaries (termed "sub-dictionaries") as values. These
sub-dictionaries have port names as keys and store names as values. For
example, the topology for the ATP example we have been considering might
look like this:

.. code-block:: python

    {
        'sodium_pump': {
            'cytoplasm': 'cell',
            'extracellularSpace': 'bloodVessel',
        },
        'metabolism': {
            'cytoplasm': 'cell',
        },
    }

-----------------
Example Composite
-----------------

To put all this information together, let's take a look at an example
composite that combines transport, growth, division, and expression
processes. This example comes from
:py:func:`vivarium.compartments.growth_division.compose_growth_division`.

.. code-block:: python

    def compose_growth_division(config):

        # declare the processes
        transport = ConvenienceKinetics(get_glc_lct_config())
        growth = Growth(config)
        division = Division(config)
        expression = MinimalExpression(config)

        # place processes in layers
        processes = [
            {'transport': transport,
             'growth': growth,
             'expression': expression},
            {'division': division}]

        # make the topology.
        topology = {
            'transport': {
                'internal': 'cell',
                'external': 'environment',
                'exchange': 'exchange',
                'fluxes': 'null',
                'global': 'global'},
            'growth': {
                'global': 'global'},
            'division': {
                'global': 'global'},
            'expression': {
                'internal': 'cell',
                'external': 'environment',
                'concentrations': 'cell_concs'}}

        # add derivers
        derivers = get_derivers(processes, topology)
        processes.extend(derivers['deriver_processes'])  # add deriver processes
        topology.update(derivers['deriver_topology'])  # add deriver topology

        # initialize the states
        states = initialize_state(processes, topology, config.get('initial_state', {}))

        options = {
            'name': 'growth_division_composite',
            'environment_port': 'environment',
            'exchange_port': 'exchange',
            'topology': topology,
            'initial_time': config.get('initial_time', 0.0),
            'divide_condition': divide_condition}

        return {
            'processes': processes,
            'states': states,
            'options': options}

You may have noticed some unfamiliar code in the above example. First,
notice that when the processes are named, they are arranged in layers.
Each layer is defined as a dictionary in a processes list. You can think
of layers as describing which processes run concurrently, and the layers
are run in order. We implement this by applying updates in-between
layers, so two processes in the same layer will operate on the same
state of the model even though one's update is actually computed first.
Processes in the second layer, though, see the model state after the
updates from the first layer have been applied.

.. WARNING:: We will soon be removing layers from Vivarium and instead
   run each process at its own timestep.

.. todo:: Update with layers removed.

Second, let's discuss derivers. Derivers let us compute information from
the model state that is useful for many processes to access. For
example, we store the mass and volume of the cell in the ``global``
store and compute it with derivers. This ``global`` store is special and
specifically for derivers. It contains information that is *computed*
from the state, but it is not directly updated by processes.

Third, we discuss the initialization of the states. This line will
appear in each composite. The inner workings of
:py:func:`vivarium.core.process.initialize_state` are beyond
the scope of this guide.

Lastly, we provide extra information in ``options``, for example the
composite name.

.. todo:: Define the available options
