Setting up your analysis
========================


Once raw data files are organized following requirements in
:doc:`data-structure`, analysis can get started.
A first step is to follow the guidelines
in :doc:`tutorial`. Here we go into more detail about:

* how to parse your data,
* how to define the observable to look at, and
* how to define conditions.

.. contents:: Contents
   :depth: 2
   :local:

Experiment and filters
------------------------

To start the analysis, you need to tell tunacell which experiment to analyse,
and whether to apply filters.


Loading an experiment
'''''''''''''''''''''

To set up the experiment, you have to give the path to the experiment folder on
your computer. We will denote this path as :code:`<path-to-exp>`, then use::

   from tunacell import Experiment
   exp = Experiment(<path-to-exp>)

By default, no filter is applied. But it is possible to associate a set
of filters to an experiment, giving instructions to how data structures will
be parsed.

.. _filter-label:

Defining the statistical ensemble
'''''''''''''''''''''''''''''''''

The statistical ensemble is the set of cells, lineages, colonies, and containers
that are parsed to compute statistics. 
In some cases, you might be likely to remove outliers, such as cells carrying
anomalous values.

To do so a :class:`FilterSet` instance must be defined
and associated it with the :class:`Experiment` object.
Detailed description about how to define filters and filter sets is in
:doc:`filters`. Here we give a simple, concrete example. Suppose you'd like to
filter out cells that didn't divide symmetrically. To do so, you first instantiate
the :class:`FilterSymmetricDivision` class::

	from tunacell.filters.cells import FilterSymmetricDivision
	myfilter = FilterSymmetricDivision(raw='length', lower_bound=0.4, upper_bound=0.6)

:code:`length` is used as the raw quantifier (assuming you have a column
:literal:`length` in your data files). Such filter requires that the
daughter cell :literal:`length` at birth must be bound within 40 and 60 percent
of the mother cell's :literal:`length` at division. Then::

	from tunacell import FilterSet
	myfset = FilterSet(filtercells=myfilter)

In the last line, the keyword argument specifies :code:`filtercells` since our
filter :code:`myfilter` acts on :class:`Cell` instances. You can define one filter
for each type of structures: :class:`Cell`, :class:`Colony`, :class:`Lineage`,
and :class:`Container`.

Once that a :class:`FilterSet` instance is defined, load it with::

	exp.set_filter(fset)

.. note::
   Filtering cell outliers may affect the tree structure, decomposing original
   tree in multiple subtrees where outlier node has been removed. Hence the
   number of trees generated from one container file depends on the filter
   applied to cells.


Defining particular samples
'''''''''''''''''''''''''''

All samples from an experiment are used for statistics, under the filtering
assumption discussed above. However, visualization of trajectories is performed
over a subset of reasonable size: this is what we'll be calling
small samples.

Small samples can be chosen specifically by user (*"I am intrigued by this
cell, let's have a look on its trajectory"*), or randomly. To do so::

   from tunacell import Parser
   parser = Parser(exp)

Note that a sample is identified by a couple of labels: the container label,
and the cell identifier. For example::

	parser.add_sample({'container_label': 'FOV_001', 'cellID': 12})

or synonymously::

	parser.add_sample(('FOV_001', 12))

This information is stored under the :attr:`samples`, and you can get a
print of the registered samples with::

	print(parser.info_samples())

You can also add randomly chosen samples::

	parser.add_sample(10)

adds 10 such samples.

Please refer to :doc:`../api/parser` for more information about how to use it.

Iterating through samples
'''''''''''''''''''''''''''''

The :code:`exp` object provides a set of iterators to parse data at each level,
with the appropriate applied filters:

* :class:`Container` level with the method :meth:`iter_containers`,
  filtered at the container level,
* :class:`Colony` level with the method :meth:`iter_colonies`,
  filtered at the container, cell, and colony levels,
* :class:`Lineage` level with the method :meth:`iter_lineages`,
  filtered at the container, cell, colony, and lineages levels
* :class:`Cell` level with the method :meth:`iter_cells`,
  filtered at the container, cell, colony, and lineages levels.

The idea behind :literal:`tunacell` is to decompose colonies into sets of lineages, *i.e.* into
sets of sequences of parentally linked cells. This way, it is possible
to extract time-series that span time ranges larger than single cell cycles.

.. note::
   Decomposition in lineages is performed randomly: at cell division,
   one daughter cell is chosen randomly to be the next step in the lineage.
   This way, lineages are independent: a given cell belongs to one, and only one
   independent lineage.

Iterating over listed samples
.............................

Use above-mentioned methods on the :class:`Parser` instance.

See :doc:`../api/parser` for more details.

Iterating over all samples
...........................

Use above-mentioned methods on the :class:`Experiment` instance.

See :doc:`../api/experiment` for more details.


Defining the observable
-----------------------

To define an observable, *i.e.* a measurable quantity that evolves through time,
use the :class:`Observable` class::

    from tunacell import Observable

and instantiate it with parameters to define a particular observable.

First parameter is the name to give to the observable (to find it back in the
analysis process).

Second, mandatory parameter is the column to use as raw data
(*e.g.* 'length', 'size', 'fluo', ...).

Then, it is possible to use time-lapse data (as stored in data files, or
processed using a time-derivative estimate) or to
determine the value of said raw observable at a particular cell cycle stage,
for example length at birth.

Indicating raw data 
'''''''''''''''''''''

First, one needs to indicate which column to be used in the raw data
file, by specifying :code:`raw='<column-name>'`.

When raw data is expected to be steady, or to be a linear function of time
within cell cycle, then use :code:`scale='linear'` (default setting). When it is
expected to be an exponential function of time within cell cycle, use
:code:`scale='log'`. We will mention below how this parameter affects some
procedures.

Raw data can be used as is, or further processed to provide user-defined observable. 
Two main modes are used to process raw data:

* The *dynamics* mode is used when one wants to analyze observables for all time
  points; examples are: length, growth rate, ...
* The *cell-cycle* mode indicates observables that are defined as a single value
  per cell cycle; examples are: length at birth, average growth rate, ...

Dynamic mode
'''''''''''''

It corresponds to the parameter :code:`mode='dynamics'`.
It sets automatically the *timing* parameter as :code:`timing='t'` where ``t``
stands for time-lapse timing. It is meant to study observables for all time
points (time-lapse, dynamic analysis).


Cell-cycle modes
''''''''''''''''

Cell-cycle modes are used when one observable need to be quantified at the cell-cycle
level, *i.e.* quantified once per cell cycle.There are few cell cycle modes:

* :code:`mode='birth'`: extrapolates values to estimate observable at cell birth;
* :code:`mode='division'`: extrapolates values to estimate observable at cell
  division;
* :code:`mode='net-increase-additive'`: returns the difference between division
  and birth values of observable;
* :code:`mode='net-increase-multiplicative'`: returns the ratio between division
  and birth values of observable;
* :code:`mode='average'`: returns the average value of observable along cell
  cycle;
* :code:`mode='rate'`: proceeds to a linear/exponential fit of observable
  depending on the chosen ``scale`` parameter. In fact, the procedure always
  performs linear fits, when :code:`scale='log'` the log of raw data is used,
  thereby performing an exponential fit on raw data.


Choosing the timing
'''''''''''''''''''

For dynamic mode, the only associated timing is ``t`` (stands for "time-lapse").
The parameter ``tref`` may be used to align time points. When provided as a
number, it will be substracted to acquisition time. A string code can be given,
``'root'`` that aligns data with the colony's root cell division time (caution:
when filtering happens, some cells that were acquired at the middle of your
experiment can become root cells if their parent cell is an outlier; this may
affect dangerously the alignement of your time-series).

For cell-cycle modes it associates
to the estimated observable a time-point to be chosen between:

* ``b``: time at birth, when known;
* ``d``: time at division, when known;
* ``m``: time at mid-point trhough cell-cycle;
* ``g``: generation index, which can be used in conjunction with the parameter
  ``tref``. When the later is set to a floating number, generation index will
  be offset to the generation index of the cell's ancestor that lived at this
  time of reference if it exists, otherwise, data from this lineage is discarded
  in analysis. When :code:`tref=None`, then the generation index is relative to
  the colony to which belongs current cell.

End-point values are estimated by extrapolation. This is because cell divisions
are recorded halfway between parent cell last frame and daughter cell first
frame. The extrapolation uses local fits over ``join_points`` points.

.. warning::

   generation index may be used with care in statistical estimates over the
   dynamics of dividing cells, since generation 0 for a given colony
   does not necessarily correspond to generation 0 of another colony.


Differentiation
'''''''''''''''

In *dynamics* mode, differentiation is obtained either by default using finite
differences with two consecutive points, either by a sliding window fit.
For an observable :math:`x(t)`, depending on the chosen scale, linear or log,
it returns the estimate of :math:`\frac{dx}{dt}` or
:math:`\frac{d}{dt} \log x(t)` respectively.


Local fit estimates
'''''''''''''''''''

As finite difference estimates of derivatives are very sensitive to
measurement precision, the user can opt for a local fitting procedure.

This procedure can be applied to estimate derivatives, or values of the
observables by performing local linear fit of the scaled observable over
a given time window. To use said option, user needs to provide the time window
extent, *e.g.* :code:`time_window=15`, will proceed to a local fit over
a time window of 15 units of time (usually minutes).

Such a local fit procedure restricted to scanning cell-cycle time segments
would lead to a loss of exploitable times, as large as the time window,
for each cell. To deal with that, the procedure provide a way to use daughter
cell information to "fill data estimates" towards the end of cell-cycle.
The parameter :code:`join_points=3` indicates that end-point values are
estimated using 3 first frames, or 3 last frames.

.. warning::

   Using local fitting procedure is likely to artificially correlate time points
   over the time window time range. Such option can help with data visualization
   since it smoothens measurement errors, but **extreme caution** is adviced when
   this feature is used in statistical analysis.

Examples
''''''''

Let's assume that raw data column names include ``'length'`` and ``'fluo'``.

Example 1: length vs. time
..........................

This is the trivial example. We stay in dynamic mode, and we do not associate
any further processing to collected data::

	>>> length = Observable(name='length', raw='length')

Example 2: length at birth
..........................

We go to the corresponding cell-cycle mode with the appropriate timing::

	>>> length_at_birth = Observable(name='birth-length', raw='length', mode='birth', timing='b')

.. note::

   one could associate the observable length at birth with another timing,
   *e.g.* time at mid cell cycle.

Example 3: Fluorescence production rate (finite differences)
.............................................................

::

	>>> k = Observable(name='prod-rate', raw='fluo', differentiate=True)

Example 4: Fluorescence production rate (local fit)
...................................................

We found that the later led to really noisy timeseries, so we choose to produce
local estimates over 3 points, in an experiment where acquisition period is
4 minutes, it means to have a 12 minutes time-window::

	>>> kw = Observable(name='window-prod-rate', raw='fluo', differentiate=True, scale='linear',
                            local_fit=True, time_window=12.)

It computes

.. math::

	\frac{d}{dt} \mathrm{fluo}(t)

using 12 minutes time windows.

Example 5: Fluorescence production rate (using local fits) at birth
...................................................................

And we want to have it as a function of generation index, setting 0 for cells
that live at time 300 minutes::

	>>> kw = Observable(name='window-prod-rate-at-birth'raw='fluo', differentiate=True, scale='linear',
                            local_fit=True, time_window=12.,
	                    mode='birth', timing='g', tref=300.)


Conditional analysis
--------------------

We saw in :ref:`filter-label` that one can define filters that act on
cells, or colonies, and to group them in a :class:`FilterSet` instance that
essentially sets the statistical ensemble over which analysis is performed.

There is another utility of these :class:`FilterSet` objects: they may define
sub-ensembles over which analysis is performed in order to compare results
over chosen sub-populations. One example is to "gate" cell-cycle quantifiers
and observe the statistics of the different sub-populations. Here we extend
the gating procedure to analyse any dynamic observable.

To do os, a list of :class:`FilterSet` instances, one per condition, can be
provided to our analysis functions. We refer to the following users pages for
further reading on how to use filters, see :doc:`filters`, and how to run
statistical analysis :doc:`statistics`.
