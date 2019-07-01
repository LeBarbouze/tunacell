Introduction to ``tunacell``
============================

``tunacell`` is a computational tool to study the dynamic properties of growing
and dividing cells. It uses raw data extracted from time-lapse microscopy and
provides the user with a set of functions to explore and visualize data, and
to compute the statistics of the dynamics.

It has been designed with the purpose of studying growth and gene expression
from fluorescent reporters in *E. coli*, but should be applicable to any
growing and dividing micro-organism.

``tunacell`` operates after segmentation and tracking has been performed on images.
For this first step, there are many tools available on the web,
see :ref:`segmentation` for a short list,
and we will try in the future to facilitate the segmentation output to
tunacell's input pipeline for some of these tools.

What it does
------------

``tunacell`` reads segmented/tracked output and is able to reconstruct lineages
and colony structures from which it can directly

1. plot the trajectories of user-defined observables for small samples,
   to gain qualitative intuition about the dynamic processes;
2. perform statistical analysis of the dynamics, to gain some quantitative
   insights about these processes.

It also provides a Python API for the user to tailor
his/her specific analysis.

How it works
------------

``tunacell`` is not as smart as you are, but with the appropriate parameters
it will compute information faster than you can.

In order to perform any analysis, the user has to define the following
quantities:

* the particular observable you want to look at,
  whether it is directly time-lapse raw values (e.g. length, total cell
  fluorescence), differentiation of it (e.g. growth rate, production rate),
  or quantities defined at cell-cycle stage (e.g. birth length,
  cell-cycle increase in fluorescence);
* the set of filters that will define the ensemble over which statistics are
  computed;
* the set of conditions that may define subgroups of samples over which
  comparison of results is relevant.

These steps are described in :doc:`users/settings`.

After these quantities are defined, API highest level functions are designed to

* plot trajectories as visual examples (see :doc:`users/plotting-samples`);
* compute the statistics of the dynamics (see :doc:`users/statistics`);
* visualize the results of such computations (idem).

Lower-level API-functions may be used by experienced users to tailor
new, specific analyses.

Why should I use ``tunacell``?
------------------------------

Because as a Python enthusiast, how cool is to say to your colleagues at
the next conference you'll attend: "I use both Python and tuna to analyze data
about how bacteria struggle in life" [#f1]_.

One of the novelties of ``tunacell`` is to provide a powerful tool to perform
conditional analysis of the dynamics. By conditional, we mean performing
statistical computations over user-defined subgroups of the original sample
ensemble.

A first set of functions to make these subgroups are already defined.
Although it is not exhaustive, the pre-defined set of subgrouping functions
already covers a wide range of possibilities.

Import/export functions have been implemented to save analyses altogether with
their parameters to allow the user to keep a structured track of what has been
done, making easier to collaborate on the analysis step.

Finally, experienced users will find it useful to be able to extend tuna's
framework, by designing new filtering functions, or implementing statistical
analyses tailored to their particular project.


Where to start then?
--------------------

We encourage readers to start with the :doc:`users/tutorial` that will present
the features described above on a simple, numerically generated dataset.

Then plug-in your data, check the documentation, and discover how cool
micro-organisms are (on a dynamical point of view).


Contribute
-----------

If you find any bug, feel free to report it.

We also welcome contributors that point directly solutions to bugs, or that
have implemented other functions useful for analysing the dynamics of growing
micro-organisms.

.. _segmentation:

Segmentation and tracking tools
-------------------------------

Before using ``tunacell`` you need to have segmented and tracked images from your
time-lapse movie. Many different softwares exist, depending on your
experimental setup. Some exemples are listed (non-exhaustively):

* `SuperSegger`_ : a Matlab-based software, with GUI. Uses machine learning
  principles to detect false divisions. Adapted for micro-colony growth in agar
  pads-like setups. Segment brightfield images (possible to inverse fluorescent
  images).
* `oufti`_ is also a Matlab-based software, following the previous
  `microbetracker`_ software developed by the same group.
* `moma`_ is a Java-based software particularly adapted to mother machine-like
  setups (and their `paper`_).
* `ieee_seg`_ is another software adapted to mother machine-like setups.
 

.. _SuperSegger: http://mtshasta.phys.washington.edu/website/SuperSegger.php
.. _oufti: http://www.oufti.org/
.. _microbetracker: http://microbetracker.org/
.. _moma: https://github.com/fjug/MoMA
.. _paper: http://biorxiv.org/content/early/2016/09/20/076224
.. _ieee_seg: http://ieeexplore.ieee.org/document/7299289/?reload=true


.. rubric:: Footnotes

.. [#f1]
   It is highly recommended to double check the conference topic beforehand.


..
    Package overview
    ----------------

    ``tunacell`` is designed to perform the following:

    *   read flat files where time-lapse data is stored
    *   re-build lineages/trees
    *   load dynamics data in memory for further analysis
	    - possibility to filter data,
	    - compute and store user-defined data,
	    - local fitting/smoothing procedures...
    *   perform statistical analysis of the dynamics:
	    - expectation values, variance estimates,
	    - auto- and cross-correlation analysis,
	    - conditional analysis...
    *   data visualization (small samples, statistics)
    *   perform numerical simulations on dividing cells

