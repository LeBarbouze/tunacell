      _                              _ _ 
     | |_ _   _ _ __   __ _  ___ ___| | |
     | __| | | | '_ \ / _` |/ __/ _ \ | |   tune-a-cell
     | |_| |_| | | | | (_| | (_|  __/ | |   Analysis of Timeseries from dividing UNicellular microorganisms
      \__|\__,_|_| |_|\__,_|\___\___|_|_|
     

`tunacell`, or in its shorter form `tuna`,
is a Python package that provides tools to analyze data from time-lapse
movies of dividing micro-organisms.

# Main features

* Analysis of time-series defined over many cell-cycles
* Computation of average values, variance, number of samples over time (easy)
* Computation of autocovariance for a single observable, and cross-covariance
  for a couple of observables over time (less easy)
* Plotting functions of the computed statistics
* Filtering-out cells
* Conditional analysis (computation of statistics among subgroups of cells,
  lineages, colonies, or containers (*aka* fields of view)
* Data visualization of small samples (looking at both the dynamic, and
  the tree structure)
* Export of computed data as text files in a comprehensive folder structure

# Install

tunacell has been developed with Python 2.7, and is now compatible with Python 3.6
(note that there might be hidden compatibility bugs).
It can be installed locally using pip.
Clone (or download) the repo and make a local (editable) install using pip:

    pip install -e .

## Local install, and `virtualenv`

If Python is installed system-wide, you may have to `sudo` the command above.
When it's not possible you may give the option to install it on the user
directory:

   pip install -e --user .

A better solution when Python, pip, and virtualenv are installed on the system,
is to create a virtual environment where you plan to work with tunacell.
The Makefile does the job, run the command:

   make virtualenv

and then

   source venv/bin/activate

to activate the virtual environment. Then you can run the `pip install -e.`
command, without worrying about permissions since everything will be installed
locally, and accesses only when your virtual environment is active.
When you finish working with tunacell, type:

   deactivate

and that's it.


## Dependencies

tunacell depends on few libraries that are automatically installed if you are
using pip.

[Numpy][], [Scipy][], [matplotlib][] are classic libraries,
as well as [pandas][] that is used
to provide the user with DataFrame objects for some statistical analyses.

The tree-like structure arising from dividing cells
has been implemented using the [treelib][] library.

[Scipy]: http://www.scipy.org/ "The Scipy package"
[Numpy]: https://docs.scipy.org/doc/numpy-dev/user/index.html "Numpy"
[pandas]: http://pandas.pydata.org/ "pandas"
[matplotlib]: http://matplotlib.org/ "matplotlib"
[treelib]: https://github.com/caesar0301/treelib  "Treelib library"

## New to Python

Python is a computer programming language and to use it, you need a Python
interpreter. To check whether you have the Python interpreter installed
on your system, run the following command in a terminal:

    python -V

If the answer shows something like ``Python 2.7.x``, or
``Python 3.6.y``, you're good to go.
Otherwise you should install it, either directly downloading
[the source files][python-downloads],
or using a friendlier package that will guide you,
such as [anaconda][anaconda].

After that you should be ready, and pip should be automatically installed. Again
try:

    pip -V

If it is not installed, you may check [this to install pip][install-pip].

Then get back to install instructions above.

[python-downloads]: https://www.python.org/ "Python"
[install-pip]: https://pip.pypa.io/en/stable/installing/ "Install pip"
[anaconda]: https://docs.continuum.io/ "Anaconda"

# Documentation

Documentation can be found [here][tunadocs], with an [introduction][tunadocs-intro],
a [quick demo tutorial][tunadocs-tutorial], a user manual that guides you
through tunacell utilization.

# Using, and learning to use tunacell

There are two main options to dive into tunacell.

First option is to go through the user manual of the documentation, stepping
through each point. Although it might give you an in-depth, logical introduction
to tuna, it might be tedious as a first approach.

The second, pragmatic option is to run **sequentially** the scripts in the ``scripts``
folder. In a terminal:

    cd <your-location-for-tunacell-repo>/scripts
    python simurun.py  # creates by default the simutest experiment
    python tutorial.py
    python plotting-samples.py
    python univariate-analysis.py
    python univariate-analysis-2.py
    python bivariate-analysis.py

Once you've run, read, and understood the bits of code in these files, consider
yourself as a tunacell expert.

If you got how it works, plug your data in
(look at [how to format input files][tunadocs-data-structure])
and use tunacell API to write your
scripts and discover new things about the dynamical properties of your cells!

# Development

## About this version

Current version is now ready for public in alpha testing. 
Bugs may come up quickly,
please report them with an Issue, or better, fork, make the patch, and PR :)

## Added features

* Made the dynamic analysis API a bit clearer (hopefully)
* Added few scripts to introduce tunacell API
* `FunctionalObservable` class has been added in tuna/base/observable.py: 
  it allows the user to define a new observable as a function of other
  observables. For instance, it can be helpful when one wants to rescale a
  dynamic, time-lapse observable (say, growth rate) by a cell-cycle observable
  (say, birth growth rate).
* Noticeable change how hidden observables are collected to perform all computation.
  A hidden observable is an observable involved in filtering (filter out those
  cells that have abnormal values for this observable), or present as an
  argument of a `FunctionalObservable`. They need to be computed internally and
  there is a distinction between regular `Observable` and the new
  `FunctionalObservable` classes. It changed the parameterization of
  `Lineage.get_timeseries` method (that can be used to parse data, see the end
  of `bivariate-analysis.py` script for a short example). 
  `Cell` instances are now equiped with protection against building mutltiple
  times a given observable (although the protection needs to be specified,
  it is not automatic).

## Future work

- [x] Make tunacell Python 3 compatible
- [ ] Add test cases (coverage is really poor right now)
- [ ] Add features for static statistical analysis (distributions, scatter-plots, ...)
- [ ] Add GUI

[tunadocs]: http://www.joachimrambeau.com/pages/_tunadocs/index.html "Tunacell documentation"
[tunadocs-intro]: http://www.joachimrambeau.com/pages/_tunadocs/intro.html "Introduction to tunacell"
[tunadocs-tutorial]: http://www.joachimrambeau.com/pages/_tunadocs/tutorial.html "10 minute tutorial"
[tunadocs-data-structure]: www.joachimrambeau.com/pages/_tunadocs/docs/_build/html/users/data-structure.html "Tunacell input format"

