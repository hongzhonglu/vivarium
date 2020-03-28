========
Glossary
========

.. glossary::

    ABM
    Agent-Based Model
    Agent-Based Models
        Agent-based modeling is a modeling paradigm in which
        population-level phenomena are emergent from interactions among
        simple agents. An agent-based model is a model constructed using
        this paradigm.

    WCM
    Whole-Cell Model
    Whole-Cell Models
        Whole-cell models seek to simulate a cell by modeling the
        molecular mechanisms that occur within it. For example, a cell's
        export of antibiotics might be modeled by the transcription of
        the appropriate genes, translation of the produced transcripts,
        and finally complexation of the translated subunits. Ideally the
        simulated phenotype is emergent from the modeled processes,
        though many such models also include assumptions that simplify
        the model.

    MSM
    Multiscale Model
    Multiscale Models
        Multiscale models use different spatial and temporal scales for
        their component sub-models. For example, Vivarium models a
        cell's internal processes and the interactions between cells and
        their environment at different temporal scales since these
        processes require different degrees of temporal precision.

    Process
    Processes
        A process in Vivarium models a cellular process by defining how
        the state of the model should change at each timepoint, given
        the current state of the model. During the simulation, each
        process is provided with the current state of the model and
        the timestep, and the process returns an update that changes
        the state of the model.

    Timestep
    Timesteps
        The amount of time elapsed between two timepoints. This is the
        amount of time for which processes compute an update. For
        example, if we discretize time into two-second intervals, then
        each process will be asked to compute an update for how the
        state changes over the next two seconds. The timestep is two
        seconds.

    Timepoint
    Timepoints
        We discretize time into timepoints and update the model state at
        each timepoint. We collect data from the model at each
        timepoint. Note that each compartment may be running with
        different timesteps depending on how finely we need to
        discretize time.

        .. todo:: How does this work with the returned timeseries data?

    Composite
    Composites
        Composites model how a set of :term:`processes` interact through
        shared :term:`stores`. When the composite is simulated, each
        store gets updated by one or more of the processes in the
        composite. This lets us simulate how the state of the model
        evolves as all the processes run simultaneously.

    Store
    Stores
        The state of the model is broken down into stores, each of which
        represents the state of some physical or conceptual subset of
        the overall state. For example, a cell model might have a store
        for the proteins in the cytoplasm, another for the transcripts
        in the cytoplasm, and one for the transcripts in the nucleus.
        Each :term:`variable` must belong to exactly one store.

    Variable
    Variables
        The state of the model is a collection of variables.  Each
        variable stores a piece of information about the full model
        state. For example, the concentration of glucose in the
        cytoplasm might be a variable, while the concentration of
        glucose-6-phosphate in the cytoplasm is another variable. The
        extracellular concentration of glucose might be a third
        variable. As these examples illustrate, variables are often
        track the amount of a molecule in a physical region. Exceptions
        exist though, for instance whether a cell is dead could also be
        a variable.

    Port
    Ports
        When a :term:`process` needs access to part of the model state,
        it will be provided a :term:`store`. The ports of a process are
        what the process calls those stores. When running a process, you
        provide a store to each of the process's ports. Think of the
        ports as physical ports into which a cable to a store can be
        plugged.

    Compartment
    Compartments
       TODO
       .. todo:: Write compartment glossary entry

    Topology
    Topologies
        A topology defines how :term:`stores` are associated to
        :term:`ports`. This tells Vivarium which store to pass to each
        port of each process during the simulation.

    Masking
        When Vivarium passes stores to processes, it includes only the
        variables the process has requested. We call this filtering
        masking.
