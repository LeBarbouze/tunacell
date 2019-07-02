Tunacell's data model
=====================

tunacell's top level data structure matches input files scaffold. Raw data is stored
in :class:`Cell` instances connected through a tree structure arising from
cell divisions.

Top level structures: :class:`Experiment` and :class:`Container`
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

tunacell's top level structure is :class:`Experiment` and handles the experiment.
We refer to the API documentation for details about attributes and methods.
In particular, it stores the list of container files that allows to open/read
such containers.

These are stored under :class:`Container` instances, which label is set by file
name. Such objects gets two major attributes: :attr:`cells` and :attr:`trees`.
The former is the list of :class:`Cell` instances imported from raw data file,
the latter is the list of reconstructed trees formed by dividing cells, stored
as :class:`Colony` instances.

Low-level structures: :class:`Cell` and :class:`Colony`
'''''''''''''''''''''''''''''''''''''''''''''''''''''''

These classes are derived from the `treelib`_ package :class:`Node` and
:class:`Tree` classes respectively.

Raw data is stored under the :attr:`data` attribute of :class:`Cell` instances.

Methods are defined at the :class:`Container` level to retrieve objects
corresponding to an identifier. More importantly there is an iterator over
colonies that can be used when parsing data for statistical analysis.

Tree decomposition: :class:`Lineage`
''''''''''''''''''''''''''''''''''''

For studying dynamics over times larger than one or few cell cycles, it is
necessary to build timeseries of observables over sequences of more than one
cells.

We use features from the treelib package to decompose trees in independent
lineages. A lineage is a sequence :math:`{c_i}_i` of cells related through
successive divisions:
cell :math:`c_i` is a daughter of cell :math`c_{i-1}`, and the mother of
cell :math:`c_{i+1}`.

One way to decompose a tree in lineages is to build the sets of lineages from
root to all leaves. Such decomposition implies that some cells may belong to
more than one lineage. Using such decomposition require some
statistical weighting procedure.

To avoid such weighting procedure, we used a decomposition in independent
lineages. Such decomposition ensures that each cell is counted once and only
once. More specifically our method to decompose a tree in independent lineages
is to traverse the tree starting from the root and choosing randomly one
daughter cell at each division until a leaf is reached, repeatidly.

A lineage is defined as a :class:`Lineage` instance. Such object gets
method to build the corresponding timeseries for a given observable.

.. _treelib: https://github.com/caesar0301/treelib

