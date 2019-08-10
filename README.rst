===========
tunacell
===========

|Build-Status| |Documentation-Status| |PyPI-Version| |License|

.. |Build-Status| image:: https://travis-ci.com/LeBarbouze/tunacell.svg?branch=develop)
   :target: https://travis-ci.com/LeBarbouze/tunacell
   :alt: Travis CI Status

.. |Documentation-Status| image:: https://readthedocs.org/projects/tunacell/badge/?version=latest
   :target: https://tunacell.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status

.. |PyPI-Version| image:: https://img.shields.io/pypi/v/tunacell?color=blue
   :target: ttps://pypi.org/project/tunacell/

.. |License| image:: https://img.shields.io/github/license/LeBarbouze/tunacell


``tunacell`` is a Python package that provides tools to analyze data from time-lapse
movies of dividing micro-organisms.

Main features
=============

- Analysis of time-series defined over many cell-cycles
- Computation of average values, variance, number of samples over sampling time or generations
  (the former is better used for dynamic quantities such as instantaneous growth rate,
  the latter for cell-cycle observables such as birth size)
- Computation of autocovariance for a single observable, and cross-covariance
  for a couple of observables over time
- Plotting functions of the computed statistics
- Filtering-out cells
- Conditional analysis: computation of statistics for subgroups of cells,
  lineages, colonies, or containers (*aka* fields of view) where the subgroups are defined
  respective to conditions (the term *gating* is employed elsewhere)
- Data visualization of small samples
- Export of computed data as text files in a comprehensive folder structure
- Numerical simulations
- Input type includes text files (comma- or tab-separated values) and output from
  segmentation software (only Supersegger's output is supported)


Install
=======

Install the package from PyPI:

    $ pip install tunacell

It is adviced to run this command in a virtual environment.


Testing the install
===================

Unit tests
----------

tunacell comes with unit tests that can be run with pytest. Run:

	make test

to execute those tests. Unit test coverage is getting better but is far from
exhaustive.

End-to-end tests
----------------

To check whether most features works appropriately, clone the
repo and run

	make full-demo

It should not raise any Python error (though some warnings may show up).

(note that these scripts writes files and folders in a new ``tmptunacell``
folder in your home directory--taking roughly 13 MB when all scripts have been
launched)


`Quick start`_
==================

.. _`quick start`: https://tunacell.readthedocs.io/en/latest/users/tutorial.html

This 10 minute tutorial gives an overview of tunacell and its API.

For a deeper, tutorial-like exploration of tunacell API,
consider running **sequentially** the scripts in
the ``scripts`` folder as it is done in the full-demo recipe in Makefile.
Just change the ``--time .5`` option to interactive mode by adding ``-i``::

   cd <your-location-for-tunacell-repo>/scripts
   $ tunasimu -s 42  # makes the simutest numerically simulated experiment
   $ python tutorial.py -i
   $ ...

(note that these scripts writes files and folders in a new ``tmptunacell``
folder in your home directory--taking roughly 13 MB when all scripts have been
launched)

Once you've run, read, and understood the bits of code in these files, consider
yourself as a tunacell expert.

If you got how it works, plug your data in
(look at [how to format input files][tunadocs-data-structure])
and use tunacell API to write your
scripts and discover new things about the dynamical properties of your cells!


Documentation_
==============

.. _documentation: https://tunacell.readthedocs.io/en/latest/

There is an extensive documentation_. Start with the `introduction to tunacell
<https://tunacell.readthedocs.io/en/latest/intro.html>`_.


Dependencies
=============

Numpy_, Scipy_, matplotlib_ are classic libraries,
as well as pandas_ that is used
to provide the user with DataFrame objects for some statistical analyses.

The tree-like structure arising from dividing cells
has been implemented using the treelib_ library.

We use pyYAML_ to parse yaml files such as metadata or other library-created
files, and tqdm_ package for progress bars.

.. _Scipy: http://www.scipy.org/
.. _Numpy: https://docs.scipy.org/doc/numpy-dev/user/index.html
.. _pandas: http://pandas.pydata.org/
.. _matplotlib: http://matplotlib.org/
.. _treelib: https://github.com/caesar0301/treelib
.. _pyYAML: https://pypi.python.org/pypi/PyYAML
.. _tqdm: https://pypi.python.org/pypi/tqdm


Contributing
============

As bugs may come up quickly, we encourage users to report them with an Issue.
We further encourage experienced users to propose patches and pull requests :)

We also welcome any suggestion for improvement (though an Issue,
or an email--see setup). Thanks for your help!


