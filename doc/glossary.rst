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
