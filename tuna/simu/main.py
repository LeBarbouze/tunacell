#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
tuna package
============

simu module
~~~~~~~~~~~

main.py module
--------------

This module collects class definitions and functions used for general
tasks when performing a numerical simulation of a process within a tree
structure.
"""
from __future__ import print_function

import numpy as np

from tuna.base.cell import CellError
from tuna.base.cell import Cell


class SimuParams(object):
    """Stores parameters for simulation parameters

    Parameters
    ----------
    nbr_container : int
        Number of simulated containers
    nbr_colony_per_container : int
        Number of colonies per container
    start : float
        time at which simulation starts. Root cell age at initial time is set
        randomly. Root cell birth time is then less or equal to start time.
    stop : float
        time at which simulation stops. Leaf cell division time is greater or
        equal to stop time.
    interval : float
        time interval at which measurements are recorded (sampling).

    Notes
    -----
    A very informal note on the number of samples:
        A sample is a tree of cells. The tree depth depends on two things:
            - the difference between start and stop times
            - the average cell cycle time (set in DivisionParams)
        The number of evaluated values for an exaclty sampled process (e.g.
        the Ornstein Uhlenbeck process) is inversely proportional to the
        period 'interval' at which data is recorded.

        The total number of operations depends:
            - exponentially on the difference 'stop - start'
            - inversely proportionally to 'interval'
            - proportionally on nbr_container
            - proportionally on nbr_colony_per_container
    """

    def __init__(self, nbr_container=1, nbr_colony_per_container=1,
                 start=0., stop=100., interval=5.):
        self.nbr_container = nbr_container
        self.nbr_colony_per_container = nbr_colony_per_container
        self.start = start
        self.stop = stop
        self.interval = interval
        # metadata helper
        self.content = [('nbr_container', nbr_container),
                        ('nbr_colony_per_container', nbr_colony_per_container),
                        ('start', start),
                        ('stop', stop),
                        ('interval', interval)]
        return

    def __str__(self):
        msg = ('Simulation parameters:\n'
               '----------------------\n')
        left_column_size = 25
        right_column_size = 10
        for key, val in self.content:
            msg += '{}'.format(key).ljust(left_column_size)
            msg += ': ' + '{}\n'.format(val).rjust(right_column_size)
        return msg


class DivisionParams(object):
    """Cell division parameters.

    Here we attribute a random value for the cell cycle duration,
    based on the two parameters `mean` and `std`. So far, the random value is
    thrown out of a gamma distribution with given mean and standard deviation.
    Division is set independently of the process simulated within each cell.

    TO IMPLEMENT: choice of distribution?

    Parameters
    ----------
    mean : float
        mean value of distribution
    std : float
        standard deviation (as sqare root of variance) of distribution
    minimum : float
        minimim value that can take the outcome
        (in practice, in order not to loose any cell in the simulation, one
        must set this value larger or equal to the period of acquisition times)
    """

    def __init__(self, mean=60., std=6., minimum=5.):
        self.minimum = minimum
        self.mean = mean
        self.std = std
        self.content = [('mean', mean),
                        ('std', std),
                        ('minimum', minimum)]
        return

    def rv(self):
        theta = self.std**2 / (self.mean - self.minimum)
        k = (self.std / theta)**2 + 1.
        val = self.minimum + np.random.gamma(k, scale=theta)
        return val

    def __str__(self):
        msg = ('Division parameters:\n'
               '--------------------\n')
        left_column_size = 5
        right_column_size = 10
        for key, val in self.content:
            msg = '{}'.format(key).ljust(left_column_size)
            msg += ': ' + '{}\n'.format(val).rjust(right_column_size)
        return msg


class EcoliError(CellError):
    pass


class EcoliDivisionError(EcoliError):
    pass


class EcoliBirthError(EcoliError):
    pass


class Ecoli(Cell):
    """Ecoli class for numerical simulations of processes on trees.

    This class is based on `Cell` class, itself being built on
    `treelib.Node` class; so it is suitable for `treelib.Tree` structures.
    Each instance require a unique identifier, a parent cell (except for root
    cells), a birth time, and a given lifetime.

    Parameters
    ----------
    identifier : str
        unique identifier
    parent : str
        label of parent instance
    birth_time : float
        time at which cell cycle starts
    lifetime :
        float, duration of cell cycle
    """

    def __init__(self, identifier=None, parent=None,
                 birth_time=None, lifetime=None):

        Cell.__init__(self, identifier)

        self._lifetime = lifetime  # duration of cell cycle
        self._division_time = None  # time at cell division

        self.parent = parent
        if parent is not None:
            self.birth_time = parent.division_time
            self.birth_value = parent.division_value
        elif birth_time is not None:
            self.birth_time = birth_time
        else:
            raise EcoliBirthError

        if self.birth_time is not None and self.lifetime is not None:
            self.division_time = self.birth_time + self.lifetime

        self._rec_times = []  # recording times

        return

    @property
    def lifetime(self):
        "Get lifetime of Ecoli instance. See below for Setter."
        return self._lifetime

    @lifetime.setter
    def lifetime(self, value):
        "Set lifetime of Ecoli instance. See above for Getter."
        self._lifetime = value
        return

    @property
    def birth_value(self):
        "Get value of simulated process when cell cycle starts."
        return self._birth_value

    @birth_value.setter
    def birth_value(self, value):
        "Set value of simulated process when cell cycle starts."
        self._birth_value = value
        return

    @property
    def division_value(self):
        "Get value of simulated process when cell cycle ends."
        return self._division_value

    @division_value.setter
    def division_value(self, value):
        "Set value of simulated process when cell cycle ends."
        self._division_value = value
        return

    @property
    def rec_times(self):
        "Get recording times. See below for setter."
        return self._rec_times

    @rec_times.setter
    def rec_times(self, ndarray):
        "Set recording times. See above for getter."
        self._rec_times = ndarray
        return
