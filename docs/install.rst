Install
=======

The easiest way to install tunacell is from wheels::

   pip install tunacell

However some introductory tutorial and scripts are missing in the library.
To get them you can visit the GitHub repository::

   https://github.com/LeBarbouze/tunacell

where you can copy/paste these scripts (look into the scripts folder).

To get everything, a good solution is to fork the repository to your local account
and/or to clone the repository on your computer. Change directory in your local
repo and do a local install::

   pip install -e .

(the :code:`-e` option stands for :literal:`editable`).
With such a clone install, the scripts are in the same place, and you can use the Makefile to run tutorials/demos.


Local install
--------------

If Python is installed system-wide, you may have to `sudo` the command above.
When it's not possible you may give the option to install it on the user
directory:

    pip install -e --user .

Virtual environment
-------------------

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


Dependencies
------------

tunacell depends on few libraries that are automatically installed if you are
using pip.

Numpy_, Scipy_, matplotlib_ are classic libraries,
as well as pandas_ that is used
to provide the user with DataFrame objects for some statistical analyses.

The tree-like structure arising from dividing cells
has been implemented using the treelib_ library.

We use pyYAML_ to parse yaml files such as metadata or other library-created
files, tqdm_ package for progress bars, and tabulate_ for fancy tabular printing.

.. _Scipy: http://www.scipy.org/ "The Scipy package"
.. _Numpy: https://docs.scipy.org/doc/numpy-dev/user/index.html "Numpy"
.. _pandas: http://pandas.pydata.org/ "pandas"
.. _matplotlib: http://matplotlib.org/ "matplotlib"
.. _treelib: https://github.com/caesar0301/treelib  "Treelib library"
.. _pyYAML: https://pypi.python.org/pypi/PyYAML "Yaml parser"
.. _tqdm: https://pypi.python.org/pypi/tqdm "tqdm progress bar"
.. _tabulate: https://pypi.python.org/pypi/tabulate "Pretty-print tabular data"

New to Python
--------------

Python is a computer programming language and to use it, you need a Python
interpreter. To check whether you have the Python interpreter installed
on your system, run the following command in a terminal:

    python -V

If the answer shows something like ``Python 2.7.x``, or
``Python 3.6.y``, you're good to go.
Otherwise you should install it, either directly downloading Python_,
or using a friendlier package that will guide you,
such as anaconda_.

After that you should be ready, and pip should be automatically installed. Again
try:

    pip -V

If it is not installed, you may check [this to install pip][install-pip].

Then get back to install instructions above.

.. _Python: https://www.python.org/ "Python"
.. _anaconda: https://docs.continuum.io/ "Anaconda"