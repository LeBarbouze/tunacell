[![Build Status](https://travis-ci.com/LeBarbouze/tunacell.svg?branch=next)](https://travis-ci.com/LeBarbouze/tunacell)

      _                              _ _ 
     | |_ _   _ _ __   __ _  ___ ___| | |
     | __| | | | '_ \ / _` |/ __/ _ \ | |   tune-a-cell
     | |_| |_| | | | | (_| | (_|  __/ | |   Analysis of Timeseries from dividing UNicellular microorganisms
      \__|\__,_|_| |_|\__,_|\___\___|_|_|
     

`tunacell` is a Python package that provides tools to analyze data from time-lapse
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

tunacell has been originally developed with Python 2.7 and is developed with
Python 3.6. Current major version 0 is compatible with both. We encourage users
to shift to Python 3 though.

The best way to install tunacell is to clone the repo and make a local,
editable install using pip:

    pip install -r requirements.txt
    pip install -e .

The first line install dependencies whereas the second line install an editable
version of the library.

This way, you will get script tutorials shipped with the library. Otherwise
the library can be installed remotely from wheels (*i.e.* without cloning the
repo):

	pip install tunacell

In this case, you should download manually the scripts from the `scripts`
folder.

## Local install, and `virtualenv`

If Python is installed system-wide, you may have to `sudo` the command above.
When it's not possible you may give the option to install it on the user
directory:

    pip install -e --user .

A better solution when Python is to create a virtual environment where you plan
to work with tunacell. It requires pip and virtualenv to be installed on your
machine. Then the Makefile does the job, run the command:

    make virtualenv

that will set up the virtual environment and install pytest and flake8 locally.
Activate the virtual environment with:

    source venv/bin/activate

Then you can run the `pip install -e.`, or `make pipinstall`
command, without worrying about permissions since everything will be installed
locally, and accesses only when your virtual environment is active.
When you finish working with tunacell, type:

    deactivate

and that's it.


## Dependencies

[Numpy][], [Scipy][], [matplotlib][] are classic libraries,
as well as [pandas][] that is used
to provide the user with DataFrame objects for some statistical analyses.

The tree-like structure arising from dividing cells
has been implemented using the [treelib][] library.

We use [pyYAML][] to parse yaml files such as metadata or other library-created
files, and [tqdm][] package for progress bars.

[Scipy]: http://www.scipy.org/ "The Scipy package"
[Numpy]: https://docs.scipy.org/doc/numpy-dev/user/index.html "Numpy"
[pandas]: http://pandas.pydata.org/ "pandas"
[matplotlib]: http://matplotlib.org/ "matplotlib"
[treelib]: https://github.com/caesar0301/treelib  "Treelib library"
[pyYAML]: https://pypi.python.org/pypi/PyYAML "Yaml parser"
[tqdm]: https://pypi.python.org/pypi/tqdm "tqdm progress bar"

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
to tunacell, it might be tedious as a first approach.

The second, pragmatic option is to run **sequentially** the scripts in
the ``scripts`` folder as it is done in the full-demo recipe in Makefile.
Just change the ``--time .5`` option to interactive mode by adding ``-i``:

    cd <your-location-for-tunacell-repo>/scripts
    tunasimu -s 42  # makes the simutest numerically simulated experiment
	python tutorial.py -i
	...

(note that these scripts writes files and folders in a new ``tmptunacell``
folder in your home directory--taking roughly 13 MB when all scripts have been
launched)

Once you've run, read, and understood the bits of code in these files, consider
yourself as a tunacell expert.

If you got how it works, plug your data in
(look at [how to format input files][tunadocs-data-structure])
and use tunacell API to write your
scripts and discover new things about the dynamical properties of your cells!


# Development

Current version is now ready for public in beta testing.

## Contributing

As bugs may come up quickly, we encourage users to report them with an Issue.
We further encourage experienced users to propose patches and pull requests :)

We also welcome any suggestion for improvement (though an Issue,
or an email--see setup). Thanks for your help!

## Testing the install

tunacell comes with unit tests that can be run with pytest. Run:

	make test

to execute those tests. Unit test coverage is getting better but is far from
exhaustive. To check whether most features works appropriately, clone the
repo and run

	make full-demo

It should not raise any Python error (though some warnings may show up).

(note that these scripts writes files and folders in a new ``tmptunacell``
folder in your home directory--taking roughly 13 MB when all scripts have been
launched)



[tunadocs]: http://www.joachimrambeau.com/pages/_tunadocs/index.html "Tunacell documentation"
[tunadocs-intro]: http://www.joachimrambeau.com/pages/_tunadocs/intro.html "Introduction to tunacell"
[tunadocs-tutorial]: http://www.joachimrambeau.com/pages/_tunadocs/tutorial.html "10 minute tutorial"
[tunadocs-data-structure]: www.joachimrambeau.com/pages/_tunadocs/docs/_build/html/users/data-structure.html "Tunacell input format"

