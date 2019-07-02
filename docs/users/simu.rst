Numerical simulations in ``tunacell``
=====================================

Using the script
----------------

When installed using pip, a ``tunasimu`` executable is provided to run
numerical simulations and save results on local directories. These
files can be used to try ``tunacell``'s analysis tools.

Such a command comes with various parameters, that will printed upon call::

    $ tunasimu -h

There is a list of (optional) parameters.


What are these simulations about?
---------------------------------

``ou`` is the value of the Ornstein-Uhlenbeck random process simulated
in each cell, ``int_ou`` is the integrated value of the random process reset
to zero at each cell birth, ``exp_int_ou`` is the exponential of the later
value. One can think of the Ornstein-Uhlenbeck as instantaneous growth rate of
the cell, and thus ``exp_int_ou`` can be associated to cell length.
