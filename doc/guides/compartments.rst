============
Compartments
============

As you consider how to model a colony of cells, you might naturally
divide the model into sub-models of the individual cells and a separate
sub-model of their environment. Further, you might perform a similar
division between the mitochondria and the cytoplasm that surrounds it.
In Vivarium, these divisions create :term:`compartments`.

.. note:: In theory you could break a model into arbitrary compartments,
    but we think it makes the most sense to have compartments represent
    spatial segregation of :term:`processes`. Thus we designed the
    framework around compartments whose processes are largely internal.
    Compartments can interact with each other, but it is more cumbersome
    to model than processes contained within a single compartment.

---------------------
Creating Compartments
---------------------

When Running a Simulation
=========================

If you worked through the :doc:`getting started guide
<../getting_started>`, you've actually already created compartments! When
you created
the environment and cell agents, those were each a compartment. If you
want to run a simulation using the framework's built-in capabilities,
you won't have to create a compartment directly.

Directly Creating Compartments
==============================

You can create compartments directly using either the
:py:class:`vivarium.core.process.Compartment` constructor or by
using a framework function.  Both are beyond the scope of this
documentation.

------------------------
Compartment Interactions
------------------------

Even though compartments represent segregated sub-models, they still
need to interact. We model these interactions using :term:`boundary
stores` between compartments. For example, the boundary store between a
cell and its environment might track the flux of metabolites between the
cell and environment compartments.

When compartments are nested, these boundary stores also exist between
the inner and the outer compartment. Thus nested compartments form a
tree whose nodes are compartments and whose edges are boundary stores. A
node's parent is its outer compartment, while its children are the
compartments within it.

Since boundary stores can also exist between compartments who share a
parent, you may find it useful to think of compartments and their
boundary stores as a bigraph (not a bipartite graph) where the tree
denotes nesting and all the edges (including those in the tree)
represent boundary stores.

-------------------
Parent Compartments
-------------------

In the processes and compartments we have discussed so far, we have
assumed that the locations of molecules tracked in :term:`stores` were
unimportant. This assumption breaks down for some parent compartments
like environments whose modeled space is too large to be homogenized by
diffusion faster than the model's :term:`timestep`. To model this
spatial heterogeneity, we employ space discretization, diffusion
processes, and multi-body physics.

.. todo:: Add code references once implementation is finalized

.. _space-discretization-lattice:

Space Discretization with Lattice
=================================

We model heterogeneous distributions of molecules throughout space by
discretizing space into a grid. For each rectangle in the grid, we track
the concentrations of all the molecules in the compartment. We also
track which rectangle each child compartment is in.

.. note:: A child compartment is always modeled as being in exactly one
    rectangle in the grid even if the cell extends into other
    rectangles. The cell's position is defined as the position of its
    midpoint.

.. note:: You may see child compartments referred to as "agents" because
    each parent compartment can be thought of as running an
    :term:`agent-based model`.

Boundary Stores in Lattice
--------------------------

We tell each cell about the concentrations of molecules using
:term:`variables` in the boundary store. When a cell imports or exports
a molecule, it stores the flux in the boundary store. The molecules are
then removed from or added to the rectangle in which the cell resides.
The flux between cells and their environment is called :term:`exchange`.

.. note:: We localize the impact of exchange on the environment to just
    the cell's immediate vicinity to allow cells to locally deplete
    resources or let extruded toxins accumulate.

Diffusion
=========

Of course, just because a cell deposits extruded molecules around itself
doesn't mean those molecules stay localized! We created processes to
model diffusion. We have two kinds of diffusion processes:

Diffusion Field
---------------

A diffusion field operates on a grid like that described above with
:ref:`lattice <space-discretization-lattice>`. The diffusion rate is
configurable. See :py:mod:`vivarium.processes.diffusion_field` for
details.

Diffusion Network
-----------------

A diffusion network models diffusion between membrane-separated regions.
The diffusion network operates on a graph whose nodes are the regions,
which are internally homogeneous, and whose edges are the membranes
through which molecules can diffuse. You can configure how quickly each
molecule can diffuse through each membrane.

In theory, a diffusion field could be modeled as a diffusion network;
however, diffusion networks are more computationally intensive to model.
Instead, diffusion networks can be used to model diffusion between a
cell and its environment through the membrane or a channel.

See :py:mod:`vivarium.processes.diffusion_network` for details.

Multi-Body Physics
==================

When cells share the same physical space, they will exclude each
other. Thermal energy from the environment also buffets the cells. We
use a multi-body physics engine to model these forces between
compartments. This process applies forces when two compartments overlap
by too much and small random forces to approximate thermal jitter.

This process is implemented in
:py:mod:`vivarium.processes.multibody_physics`.
