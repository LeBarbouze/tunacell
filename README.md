      _                              _ _ 
     | |_ _   _ _ __   __ _  ___ ___| | |
     | __| | | | '_ \ / _` |/ __/ _ \ | |   tune-a-cell
     | |_| |_| | | | | (_| | (_|  __/ | |   Analysis of Timeseries from dividing UNicellular microorganisms
      \__|\__,_|_| |_|\__,_|\___\___|_|_|
     

`tunacell`, or in its shorter form `tuna`,
is a Python package that provides tools to analyze data from time-lapse
movies of dividing micro-organisms.

`tuna`'s main functions are to: parse data, visualize small sample realizations,
and perform statistical analysis of the dynamics. Although not restricted to, it
is specially useful for the analysis of time-series spanning multiple cell
cycles, where data must be extracted from lineages, and not just single cells
separately. It deals with dynamic data (i.e. defined for each acquisition time)
and with cell-cycle level data. Important features are the computation of
autocovariance, and cross-covariance as functions of time.

# Main features

* Analysis of time-series defined over many cell-cycles
* Computation of average values, variance, number of samples over time (easy)
* Computation of autocovariance for a single observable, and cross-covariance
  for a couple of observables over time (less easy)
* Filtering-out of cells
* Conditional analysis
* Data visualization of small samples
* Export of computed data as text files in a comprehensive folder structure

# Install

tunacell has been developed with Python 2 and can be installed locally using pip.
Clone (or download) the repo and make a local (editable) install using pip:

    pip install -e .

## New to Python

Python is a computer programming language and to use it, you need a Python
interpreter. To check whether you have the Python interpreter installed
on your system, run the following command in a terminal:

    python -V

If the answer shows something like ``Python 2.7.x``, you're good to go.
Otherwise you should install it, either directly downloading
[python-downloads][the source files],
or using a friendlier package that will guide you,
such as [anaconda][anaconda].

After that you should be ready, and pip should be automatically installed. Again
try:

    pip -V

If it is not installed, you may check [install-pip][this to install pip].

[install-pip]: https://pip.pypa.io/en/stable/installing/ "Install pip"
[anaconda]: https://docs.continuum.io/ "Anaconda"

# Documentation

Documentation can be found [here][tunadocs], with an [introduction][tunadocs-intro],
a [quick demo tutorial][tunadocs-tutorial], a user manual that guides you
through tunacell utilization.

[tunadocs]: http://www.joachimrambeau.com/_tunadocs/index.html "Tunacell documentation"
[tunadocs-intro] : http://www.joachimrambeau.com/_tunadocs/intro.html "Introduction to tunacell"

# Using, and learning to use tunacell

There are two main options to dive into tunacell.

First option is to go through the user manual of the documentation, stepping
through each point. Although it might give you an in-depth, logical introduction
to tuna, it might be tedious as a first approach.

The second option is to open, read, and run the script files stored under the
``scripts`` folder in the repo. First start with the ``simurun.py`` script
from the command line to generate data:

    python simurun.py

Then have a look at ``tutorial.py`` for a first glance at tunacell scope.
Dive in with the ``univariate-analysis.py``, ``univariate-analysis-2.py``,
and ``bivariate-analysis.py`` to become an expert.

If you got how it works, plug your data in and use tunacell API to write your
scripts and discover new things about the dynamical properties of your cells!

# Development

## About this version

This version 0.0.7 is now ready for public, alpha testing. 
Bugs may come up quickly,
please report them with an Issue, or better, fork, make the patch, and PR :)

## Added features

Amongst a global reorganization, noticeable added feature is FunctionalObservable
that allows the user to define a new observable as a function of other
observables. For instance, it can be helpful when one wants to rescale a
dynamic, time-lapse observable (say, growth rate) by a cell-cycle observable
(say, brith growth rate).

## Future work

[] Make tunacell Python 3 compatible
[] Add features for static statistical analysis (distributions, scatter-plots, ...)
[] Add GUI

