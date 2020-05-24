========
Glossary
========

.. note:: All fully-capitalized words and phrases have the meanings
    specified in :rfc:`2119`.

.. glossary::

    ABM
    Agent-Based Model
    Agent-Based Models
        Agent-based modeling is a modeling paradigm in which
        population-level phenomena are emergent from interactions among
        simple agents. An agent-based model is a model constructed using
        this paradigm.
    
    Boundary Store
    Boundary Stores
        :term:`Compartments` interact through boundary stores that
        represent how the compartments affect each other. For example,
        between an environment compartment and a cell compartment, there
        might be a boundary store to track the flux of metabolites from
        the cell to the environment and vice versa.

    Compartment
    Compartments
        We organize our models into compartments, each of which is like
        an agent in an :term:`agent-based model`. Each compartment
        stores a :term:`composite` of :term:`processes` and
        :term:`stores`. Compartments can be nested and interact with
        neighbor, parent, and child compartments through :term:`boundary
        stores`. Thus, a model might contain a compartment for the
        environment that contains two child compartments for the two
        cells in the environment. For more details, see our :doc:`guide
        to compartments </guides/compartments>`.

    Composite
    Composites
        Composites model how a set of :term:`processes` interact through
        shared :term:`stores`. When the composite is simulated, each
        store gets updated by one or more of the processes in the
        composite. This lets us simulate how the state of the model
        evolves as all the processes run simultaneously. For more
        information, see our :doc:`guide to composites
        </guides/composites>`.

    Deriver
    Derivers
        Derivers run after all processes have run for a
        :term:`timepoint` and compute values from the state of the
        model. These computed values are generally stored in the
        ``global`` :term:`store`. For example, one common deriver uses
        the cell's mass and density to compute the volume.

    Emitter
    Emitters
        While a simulation is running, the current state is stored in
        :term:`stores`, but this information is overwritten at each
        timestep with an updated state. When we want to save off
        variable values for later analysis, we send these data to
        one of our emitters, each of which formats the data for a
        storage medium, for example a database or a Kafka message. We
        then query the emitter to get the formatted data.

    Exchange
    Exchanges
        The flux between a cell and its environment. This is stored in a
        :term:`boundary store`.

    Masking
        When Vivarium passes stores to processes, it includes only the
        variables the process has requested. We call this filtering
        masking.

    MSM
    Multiscale Model
    Multiscale Models
        Multiscale models use different spatial and temporal scales for
        their component sub-models. For example, Vivarium models a
        cell's internal processes and the interactions between cells and
        their environment at different temporal scales since these
        processes require different degrees of temporal precision.

    Port
    Ports
        When a :term:`process` needs access to part of the model state,
        it will be provided a :term:`store`. The ports of a process are
        what the process calls those stores. When running a process, you
        provide a store to each of the process's ports. Think of the
        ports as physical ports into which a cable to a store can be
        plugged.

    Process
    Processes
        A process in Vivarium models a cellular process by defining how
        the state of the model should change at each timepoint, given
        the current state of the model. During the simulation, each
        process is provided with the current state of the model and
        the timestep, and the process returns an update that changes
        the state of the model. Each process is an instance of a
        :term:`process class`.

        To learn how to write a process, check out
        :doc:`our process-writing tutorial</tutorials/write_process>`.
        For a detailed guide to processes, see :doc:`our guide to
        processes </guides/processes>`.

    Process Class
    Process Classes
        A process class is a Python class that defines a process's
        model. These classes can be instantiated, and optionally
        configured, to create :term:`processes`. Each process class must
        subclass either :py:class:`vivarium.core.process.Process`
        or another process class.

    Store
    Stores
        The state of the model is broken down into stores, each of which
        represents the state of some physical or conceptual subset of
        the overall state. For example, a cell model might have a store
        for the proteins in the cytoplasm, another for the transcripts
        in the cytoplasm, and one for the transcripts in the nucleus.
        Each :term:`variable` must belong to exactly one store.

    Template
    Templates
        A template describes a genetic element, its binding site, and
        the available downstream termination sites on genetic material.
        A chromosome has operons as its templates which include sites
        for RNA binding and release. An mRNA transcript also has
        templates which describe where a ribosome can bind and will
        subsequently release the transcript. Templates are defined in
        :term:`template specifications`.

    Template Specification
    Template Specifications
        Template specifications define :term:`templates` as
        :py:class:`dict` objects with the following keys:

        * **id** (:py:class:`str`): The template name. You SHOULD use
          the name of the associated operon or transcript.
        * **position** (:py:class:`int`): The index in the genetic
          sequence of the start of the genetic element being described.
          In a chromosome, for example, this would denote the start of
          the modeled operon's promoter. On mRNA transcripts (where we
          are describing how ribosomes bind), this SHOULD be set to
          ``0``.

          .. todo:: Is position 0 or 1 indexed?

        * **direction** (:py:class:`int`): ``1`` if the template should
          be read in the forward direction, ``-1`` to proceed in the
          reverse direction.  For mRNA transcripts, this SHOULD be ``1``.
        * **sites** (:py:class:`list`): A list of binding sites. Each
          binding site is specified as a :py:class:`dict` with the
          following keys:

            * **position** (:py:class:`int`): The offset in the sequence
              from the template *position* to the start of the binding
              site.  This value is not currently used and MAY be set to
              0.
            * **length** (:py:class:`int`): The length, in base-pairs,
              of the binding site. This value is not currently used and
              MAY be set to 0.
            * **thresholds** (:py:class:`list`): A list of tuples, each
              of which has a factor name as the first element and a
              concentration threshold as the second. When the
              concentration of the factor exceeds the threshold, the
              site will bind the factor. For example, in an operon the
              factor would be a transcription factor.

        * **terminators** (:py:class:`list`): A list of terminators,
          which halt reading of the template. As such, which genes are
          encoded on a template depends on which terminator halts
          transcription or translation. Each terminator is specified as
          a :py:class:`dict` with the following keys:

            * **position** (:py:class:`int`): The index in the genetic
              sequence of the terminator. 
            * **strength** (:py:class:`int`): The relative strength of
              the terminator. For example, if there remain two
              terminators ahead of RNA polymerase, the first of strength
              3 and the second of strength 1, then there is a 75% chance
              that the polymerase will stop at the first terminator. If
              the polymerase does not stop, it is guaranteed to stop at
              the second terminator.
            * **products** (:py:class:`list`): A list of the genes that
              will be transcribed or translated should
              transcription/translation halt at this terminator.
        
    Timepoint
    Timepoints
        We discretize time into timepoints and update the model state at
        each timepoint. We collect data from the model at each
        timepoint. Note that each compartment may be running with
        different timesteps depending on how finely we need to
        discretize time.

        .. todo:: How does this work with the returned timeseries data?

    Timestep
    Timesteps
        The amount of time elapsed between two timepoints. This is the
        amount of time for which processes compute an update. For
        example, if we discretize time into two-second intervals, then
        each process will be asked to compute an update for how the
        state changes over the next two seconds. The timestep is two
        seconds.

    Topology
    Topologies
        A topology defines how :term:`stores` are associated to
        :term:`ports`. This tells Vivarium which store to pass to each
        port of each process during the simulation.

    Update
    Updates
        An update describes how the model state should change due to the
        influence of a :term:`process` over some period of time (usually
        a :term:`timestep`).

    Updater
    Updaters
        An updater describes how an update should be applied to the
        model state to produce the updated state. For example, the
        update could be added to the old value or replace it.

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
