#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
tunacell package
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

from tunacell.base.cell import CellError
from tunacell.base.cell import Cell


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
    """Sets parameters for computing interdivision times.
    
    This settings is based on the linear results explained elsewhere.
    
    Parameters
    ----------
    lambda_ : float between 0 and 1
        lambda parameter (note trailing): mixing parameter between timer (0)
        and size (1). Adder is obtained with lambda_ 0.5.
    size_target : float
        X^* set the cutoff size for division. When lambda_ is 1, and there are
        no fluctuations, division is triggered when size reaches exactly X^*
        (fixed final size).
    use_growth_rate : str {'parameter', 'birth'}
        whether to use parameter growth rate, or local growth rate
    mode : str {'gamma', 'none'}
        'fixed': only the linear response is taken;
        'gamma': adds fluctuations over linear response (needs kwargs)
    relative_fluctuations : float (default 0.1)
        sets the standard deviation relative to mean value (default 10%)
    """
    
    def __init__(self, lambda_=0, size_target=2., use_growth_rate='parameter',
                 mode='gamma', relative_fluctuations=.1):
        self.lambda_ = lambda_
        self.size_target = size_target
        self.use_growth_rate = use_growth_rate
        self.mode = mode
        self.relative_fluctuations = relative_fluctuations
        self.content = [('lambda', lambda_),
                        ('size_target', size_target),
                        ('use_growth_rate', use_growth_rate),
                        ('mode', mode),
                        ('sd_to_mean', relative_fluctuations)]
        return
    
    def rv(self, birth_size=1., growth_rate=1./60*np.log(2.)):
        """Returns a random variate for interdivision time given X_0 and alpha.

        Parameters
        ----------
        birth_size : float (default 1.)
            birth size
        growth_rate : float (default 1./60*np.log(2.))
            default value corresponds to a size doubling time of one hour;
            must be used in 'per minutes'
        
        Returns
        -------
        tau : float
            interdivision time computed according to linear response with
            coefficient lambda plus fluctuations when mode is 'gamma'
        """
        tau = ((1. - self.lambda_) * np.log(2.) +
               self.lambda_ * np.log(self.size_target/birth_size)) / growth_rate
        if self.mode == 'gamma':
            tau = _gamma(mean=tau, std=self.relative_fluctuations*tau)
        return tau

    def __str__(self):
        msg = ('Division parameters:\n'
               '--------------------\n')
        left_column_size = 5
        right_column_size = 10
        for key, val in self.content:
            msg = '{}'.format(key).ljust(left_column_size)
            msg += ': ' + '{}\n'.format(val).rjust(right_column_size)
        return msg
    
def _gamma(mean=60, std=6.):
    theta = std**2 / mean
    k = (std / theta)**2 + 1.
    val = np.random.gamma(k, scale=theta)
    return val


#class DivisionParams(object):
#    """Cell division parameters.
#
#    Here we attribute a random value for the cell cycle duration,
#    based on the two parameters `mean` and `std`. So far, the random value is
#    thrown out of a gamma distribution with given mean and standard deviation.
#    Division is set independently of the process simulated within each cell.
#
#    TO IMPLEMENT: choice of distribution?
#
#    Parameters
#    ----------
#    mean : float
#        mean value of distribution
#    std : float
#        standard deviation (as sqare root of variance) of distribution
#    minimum : float
#        minimim value that can take the outcome
#        (in practice, in order not to loose any cell in the simulation, one
#        must set this value larger or equal to the period of acquisition times)
#    """
#
#    def __init__(self, mean=60., std=6., minimum=5.):
#        self.minimum = minimum
#        self.mean = mean
#        self.std = std
#        self.content = [('mean', mean),
#                        ('std', std),
#                        ('minimum', minimum)]
#        return
#
#    def rv(self):
#        theta = self.std**2 / (self.mean - self.minimum)
#        k = (self.std / theta)**2 + 1.
#        val = self.minimum + np.random.gamma(k, scale=theta)
#        return val
#
#    def __str__(self):
#        msg = ('Division parameters:\n'
#               '--------------------\n')
#        left_column_size = 5
#        right_column_size = 10
#        for key, val in self.content:
#            msg = '{}'.format(key).ljust(left_column_size)
#            msg += ': ' + '{}\n'.format(val).rjust(right_column_size)
#        return msg


class SampleInitialSize(object):
    """Initialize cell size given a target value for division and fluctuations.

    Parameters
    ----------
    size_cutoff : float (default 2.)
        sets the target division size X^* (it sets the scale)
    mode : str {'fixed', 'lognormal'}
        initial distribution of birth size; fixed corresponds to a single value
        which is half the cutoff, lognormal returns a lognormal variate with
        parameter sigma
    sigma : float (default 2.*np.log(2))
        sets the normal fluctuations of Y = log(X/ (X^*/2))
    """
    
    def __init__(self, size_cutoff=2., mode='fixed', sigma=2.*np.log(2.)):
        self.size_cutoff = size_cutoff
        self.mode = mode
        self.sigma = sigma
        self.content = [('size_cutoff', size_cutoff),
                        ('size_mode', mode),
                        ('size_sigma', sigma)]
        return

    def rv(self):
        if self.mode == 'fixed':
            return self.size_cutoff / 2.
        elif self.mode == 'lognormal':
            reduced = np.random.lognormal(mean=0., sigma=self.sigma)
            return reduced * self.division_size/2.


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
