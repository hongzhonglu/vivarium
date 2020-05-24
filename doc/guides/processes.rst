=========
Processes
=========

You should interpret words and phrases that appear fully capitalized in
this document as described in :rfc:`2119`. Here is a brief summary of
the RFC:

* "MUST" indicates absolute requirements. Vivarium may not work
  correctly if you don't follow these.
* "SHOULD" indicates strong suggestions. You might have a valid reason
  for deviating from them, but be careful that you understand the
  ramifications.
* "MAY" indicates truly optional features that you can include or
  exclude as you wish.

Models in Vivarium are built by combining :term:`processes`, each of
which models a mechanism in the cell. These processes can be combined in
a :term:`composite` to build more complicated models. Process models are
defined in a class that inherit from
:py:class:`vivarium.core.process.Process`, and these
:term:`process classes` can be instantiated to create individual
processes.  During instantiation, the process class may accept
configuration options.

.. note:: Processes are the foundational building blocks of models in
   Vivarium, and they should be as simple to define and compose as
   possible.

---------------
Process Classes
---------------

Each process class MUST implement the API that we describe below.

Class Variables
===============

Each process class SHOULD define default configurations in a
``defaults`` class variable. These defaults SHOULD then be read by the
constructor. For example:

.. code-block:: python
    
    class MyProcess:
        defaults = {
            'growth_rate': 0.0006,
        }  


Constructor
===========

The constructor of a process class MUST accept as its first positional
argument an optional dictionary of configurations. If the process class
is configurable, it SHOULD accept configuration options through this
dictionary.

In the constructor, the process class MUST call its superclass
constructor with its :term:`ports` and a dictionary of parameters. 

.. _constructor-define-ports:

Defining Ports
--------------

Ports MUST be specified as a dictionary with port names as keys and
lists of :term:`variable` names as values. These port names may be
chosen arbitrarily. Variable names are also at the discretion of the
process class author, but note that if two processes are to be combined
in a :term:`composite` and share variables through a shared
:term:`store`, the processes MUST use the same variable names for the
shared variables.

.. note:: Variables always have the same name, no matter which process
    is interacting with them. This is unlike stores, which can take on
    different port names with each process.

Passing Parameters to Superclass Constructor
--------------------------------------------

The dictionary of parameters SHOULD include any configuration options
not used by the process class. Any information needed by the process
class MAY also be included in these parameters. Once the object has
been instantiated, these parameters are available as
``self.parameters``, where they have been stored by the
:py:class:`vivarium.core.process.Process` constructor.

Example Constructor
-------------------

Let's examine an example constructor from a growth process class.

.. code-block:: python

    def __init__(self, initial_parameters={}):
        ports = {
            'global': ['mass', 'volume']}

        parameters = {'growth_rate': self.defaults['growth_rate']}
        parameters.update(initial_parameters)
        super(Growth, self).__init__(ports, parameters)

In this constructor, only one port, ``global``, is defined, from which
the process will only need the ``mass`` and ``volume`` variables. While
the default growth rate is ``0.0006``, this can be overridden by
including a ``growth_rate`` key in the configuration dictionary passed
to ``initial_parameters``.

.. note:: ``global`` is a special port used by :term:`derivers`. It
    stores information about the total model state that, like ``mass``
    doesn't fit into any store.

Default Settings
================

The process class MUST implement a ``default_settings`` method that can
be called with no arguments. This method MUST return a dictionary with
the ``state`` key for the default state. The dictionary MAY also contain
the following keys: ``emitter_keys`` for the emitter keys, ``schema``
for the schema, and ``deriver_setting`` for the deriver settings. We
describe each of these in turn:

.. _constructor-default-state:

Default State
-------------

The process class MUST provide a default value for each variable
included in its ports declaration in the constructor, with the exception
that variables whose values will be computed by :term:`derivers` do not
need a default value. These default values MUST be specified as a
dictionary whose keys are port names and whose values are dictionaries,
termed sub-dictionaries. Each sub-dictionary has keys of variable names
and values of variable values. For example, the growth process class we
have been discussing might have a default state like this:

.. code-block:: python

    {
        'global': {
            'mass': 1339  # Mass in fg
        }
    }

Here we exclude the ``volume`` variable, which is computed by a deriver.

Emitter Keys
------------

As the simulation runs, the total model state is recorded in the stores,
but this state is overwritten each timestep with an updated state. To
save data for analysis, we send variable values to an :term:`emitter`,
for example a Kafka emitter or one for a database. The emitter keys
specify which variables' values are sent to emitters for recording.
Emitter keys MUST be specified as a dictionary of the same form as the
:ref:`ports declaration dictionary <constructor-define-ports>`, but with
only the variables to be emitted.

.. _constructor-schema:

Schema
------

.. todo:: What else does the schema do?

In the schema, we define how this process class will specify
:term:`updates` for each variable. The available updaters are as
follows:

* ``accumulate`` is the default, and it specifies that the value of the
  variable in the update be added to the variable's current value when
  the update is applied.
* ``set`` specifies that the update value overwrite the current value.

The schema MUST take the form of a dictionary like the default state
dictionary, only the variable values are replaced with dictionaries that
MAY include the ``updater`` key with a value equal to the name of the
desired updater. Variables MAY be omitted, in which case they will take
on the default updater of ``accumulate``.

Each value in the schema MAY also specify a mass using the ``mass`` key.
If you are using the mass :term:`deriver`, each variable accessed by the
deriver MUST specify a mass.

Deriver Setting
---------------

:term:`Derivers` calculate metrics or perform conversions for us over
the course of the simulation, but they do not encode mechanism. For
example, we use them to calculate a cell's mass or convert between
counts and concentrations. We configure derivers with a list of
dictionaries, one dictionary for each deriver. For example:

.. code-block:: python

    deriver_setting = [
        {   # Configuration for one deriver
            'type': ...
        },
        {   # Configuration for another deriver
            'type': ...
        },
    ]

Example Default Settings
------------------------

Let's take a look at a potential ``default_settings`` method for our
growth process:

.. code-block:: python

    def default_settings(self):

        # default state
        default_state = {
            'global': {
                'mass': 1339
            }
        }

        # default emitter keys
        default_emitter_keys = {'global': ['mass']}

        # schema
        schema = {
            'global': {
                'mass': {
                    'updater': 'set'}}}

        # We can omit the deriver_setting key so long as we aren't using
        # derivers
        default_settings = {
            'state': default_state,
            'emitter_keys': default_emitter_keys,
            'schema': schema}

        return default_settings

Here, we set the mass to a default of 1339. We also choose to emit the
``mass`` variable's values and to overwrite the mass variable on update.

Next Updates
============

Each process class MUST implement a ``next_update`` method that accepts
two positional arguments: the :term:`timestep` and the current state of
the model. The timestep describes, in units of seconds, the length of
time for which the update should be computed.

State Format
------------

The ``next_update`` method MUST accept the model state as a dictionary
of the same form as the :ref:`default state dictionary
<constructor-default-state>`.

.. note:: In the code, you may see the model state referred to as
    ``states``. This is left over from when stores were called states,
    and so the model state was a collection of these states. As you may
    already notice, this naming was confusing, which is why we now use
    the name "stores."

Because of :term:`masking`, each
port will contain only the variables specified in the
:ref:`constructor's ports declaration <constructor-define-ports>`, even
if the linked store contains more variables.

.. WARNING:: The ``next_update`` method MUST NOT modify the states it is
    passed in any way. The state's variables are not copied before they
    are passed to ``next_update``, so changes to any objects in the
    state will affect the model state before the update is applied.

Update Format
-------------

``next_update`` MUST return a single dictionary, the update that
describes how the modeled mechanism would change the model state over
the specified time. The update dictionary MUST be of the same form as the
:ref:`default state dictionary <constructor-default-state>`, though
variables that do not need to be updated can be excluded.

Example Next Update Method
--------------------------

Here is an example ``next_update`` method for our growth process:

.. code-block:: python

    def next_update(self, timestep, states):
        mass = states['global']['mass']
        new_mass = mass * np.exp(self.parameters['growth_rate'] * timestep)
        return {'global': {'mass': new_mass}}

Recall from :ref:`our example schema <constructor-schema>` that we use
the ``set`` updater for the ``mass`` variable. Thus, we compute the new
mass of the cell and include it in our update. Notice that we access the
growth rate specified in the constructor by using the
``self.parameters`` attribute.

.. note:: Notice that this function works regardless of what timestep we
    use. This is important because different composites may need
    different timesteps based on what they are modeling.

Process Class Examples
======================

Many of our process classes have examples in the form of test functions
at the bottom. These are great resources if you are trying to figure out
how to use a process.

If you are writing your own process, please include these examples!
Also, executing the process class Python file should execute one of
these examples and save the output as demonstrated in
:py:mod:`vivarium.processes.convenience_kinetics`. Lastly, any top-level
functions you include that are prefixed with ``test_`` will be executed
by ``pytest``. Please add these tests to help future developers make
sure they haven't broken your process!

---------------------
Using Process Objects
---------------------

Your use of process objects will likely be limited to instantiating them
and passing them to other functions in Vivarium that handle running the
simulation. Still, you may find that in some instances, using process
objects directly is helpful. For example, for simple processes, the
clearest way to write a test may be to run your own simulation loop.

Simulating a process can be sketched by the following pseudocode:

.. code-block:: python

    # Create the process
    configuration = {...}
    process = ProcessClass(configuration)

    # Get the initial state from the process's defaults
    # This means the stores and ports are the same
    state = process.default_settings()['state']

    # Run the simulation in a loop for 10 seconds
    time = 0
    while time < 10:
        # We are using a timestep of 1 second
        update = process.next_update(1, state)
        # This is a simplified way to apply the update that assumes all
        # all variables are numbers and all updaters are "accumulate"
        for port in update:
            for variable_name, value in port.items():
                state[port][variable_name] += value
    # Now that the loop is finished, the predicted state after 10
    # seconds is in "state"

The above pseudocode is simplified, and for all but the most simple
processes you will be better off using Vivarium's built-in simulation
capabilities. We hope though that this helps you understand how
processes are simulated and the purpose of the API we defined.
