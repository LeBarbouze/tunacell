Filters
=======

Outliers may have escaped segmentation/tracking quality control tools and thus
there might be a need to apply further filtering when analysing their output
data, as tunacell does. For example, filamentous cells may have been reported in
data, but one might exclude them from a given analysis.

tunacell provides a set of user-customisable filters that allows user to define
properly the statistical ensemble of samples over which its analysis will be
performed.

In addition to removing outliers, filters are also used for conditional
analysis, as they allow to divide the statistical ensemble in sub-populations
of samples that verify certain rules.

Series of filters have already been defined for each of the following types:
cell, lineage, colony, and container. In addition boolean operations
AND, OR, NOT can be used within each type. Then filters of different types are
combined in :class:`FilterSet` instances: one is used to define the statistical
ensemble (remove outliers), and optionnally, others may be used to create
sub-populations of samples for comparative analyses.


.. contents:: Contents
   :depth: 2
   :local:


How individual filters work
---------------------------

Filters are instances of the :class:`FilterGeneral` class.
A given filter class is instantiated (possibly) with parameters, that define
how the filter work.
Then the instantiated object is callable on the object to be filtered.
It returns either :code:`True` (the object is valid) or :code:`False`
(the object is rejected).

Four main subclasses are derived from :class:`FilterGeneral`, one for each
structure that tuna recognizes: :class:`FilterCell` for :class:`Cell` objects,
:class:`FilterTree` for :class:`Colony` objects, :class:`FilterLineage` for
:class:`Lineage` objects, :class:`FilterContainer` for :class:`Container`
objects.


Example: testing the parity of cell identifier
''''''''''''''''''''''''''''''''''''''''''''''

The filter :class:`FilterCellIdparity` has been designed for illustration:
it tests whether the cell identifier is even (or odd).

First we set the filter by instantiating its class with appropriate parameter::

    >>> from tunacell.filters.cells import FilterCellIDparity
    >>> filter_even = FilterCellIDparity(parity='even')

For this filter class, there is only one keyword parameter, :code:`parity`,
which we have set to :code:`'even'`: accept cells with even identifier, rejects
cells with odd identifier.

First, we can print the string representation::

    >>> print(str(filter_even))
    CELL, Cell identifier is even

The first uppercase word in the message reminds the type of objects the filter
is acting upon. Then the message is a label that has been defined in the class
definition).

We set two :class:`Cell` instances, one with even identifier, and one odd::

    >>> from tunacell.base.cell import Cell
    >>> mygoodcell = Cell(identifier=12)
    >>> mybadcell = Cell(identifier=47)

Then we can perform the test over both objects::

    >>> print(filter_even(mygoodcell))
    True
    >>> print(filter_even(mybadcell))
    False

We also mention another feature implemented in the representation of such
filters::

    >>> print(repr(filter_even))
    FilterCellIDparity(parity='even', )

Such representation is the string one would type to re-instantiate the filter.
This representation is used by tuna when a data analysis is exported to text
files. Indeed, when tuna reads back this exported files, it is able to load the
objects defined in the exported session.
Hence, no need of remembering the precise parameters adopted
on a particular analysis: if it's exported, it can be loaded later on.

Creating a new filter
'''''''''''''''''''''

Few filters are already defined in the following modules:

* :mod:`tunacell.filters.cells` for filters acting on cells,
* :mod:`tunacell.filters.lineages` for filters acting on lineages,
* :mod:`tunacell.filters.trees` for filters acting on colonies,
* :mod:`tunacell.filters.containers` for filters acting on containers.

Within each type, filters can be combined with boolean operations (see below),
that allows user to explore a range of filters.
However a user may need to define its own filter(s), and he/she is encouraged
to do so following the general guidelines:

* define a :attr:`label` attribute (human-readable message, which was
  :code:`'Cell identifier is even'` in our previous example),
* define the :meth:`func` method that performs the boolean testing.

From the module :mod:`tunacell.filters.cells` we copied below the class definition
of the filter used in our previous example::

    class FilterCellIDparity(FilterCell):
        """Test whether identifier is odd or even"""

        def __init__(self, parity='even'):
            self.parity = parity
            self.label = 'Cell identifier is {}'.format(parity)
            return

        def func(self, cell):
            # test if even
            try:
                even = int(cell.identifier) % 2 == 0
                if self.parity == 'even':
                    return even
                elif self.parity == 'odd':
                    return not even
                else:
                    raise ValueError("Parity must be 'even' or 'odd'")
            except ValueError as ve:
                print(ve)
                return False

Although this filter may well be useless in actual analyses, it shows how to
define a filter class. Also have a look at filters defined in the
above-mentioned modules.

How to combine individual filters together with boolean operations
------------------------------------------------------------------

Filters already implemented are "atomic" filters, *i.e.* they perform one
testing operation. It is possible to combine many atomic filters of the
same type (*type* refers to the object type on which filter is applied:
cell, lineage, colony, container) by using Boolean filter types.

There are 3 of them, defined in :mod:`tuna.filters.main`: :class:`FilterAND`,
:class:`FilterOR`, :class:`FilterNOT`. The first two accepts any number of
filters, that are combined with the AND/OR logic respectively; the third accepts
one filter as argument.

With these boolean operations, complex combinations of atomic filters can be
created.


How to define a :class:`FilterSet` instance
-------------------------------------------

So far we saw how to use filters for each type of structures, independently:
cell, lineage, colony, and container.

The :class:`FilterSet` registers filters to be applied on each of these types.
It is used to define the statistical ensemble of valid samples, or to define
a condition (rules to define a sub-population from the statistical ensemble).

Explicitely, if we would like to use our :obj:`filter_even` from our
example above as the only filter to make the statistical ensemble, we would
define::

    from tunacell.filters.main import FilterSet
    fset = FilterSet(filtercell=filter_even)

(the other keyword parameters are ``filterlineage``, ``filtertree``, and
``filtercontainer``)
