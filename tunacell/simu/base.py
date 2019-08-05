#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Base classes for numerical simulations of growing and dividing cells

This module defines base classes to parameterize numerical simulations regarding:
- global simulation parameters such as duration of experiment or the number of samples,
- division control parameters (how cell division is triggered as a function of other cell observables),
- initial cell size sampling

See also
---------

:mod:`tunacell.simu.ou` that defines classes and functions to simulate a particular
growth process using the Ornstein-Uhlenbeck model.
"""
from __future__ import print_function

import numpy as np

from tunacell.base.cell import CellError
from tunacell.base.cell import Cell


class SimuParams(object):
    """Stores global parameters for simulation parameters

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
    period : float
        time interval at which measurements are recorded (sampling).
    seed : int (default None)
        seed for random number generation (given to numpy.random.seed())

    Notes
    -----
    A very informal note on the number of samples:
        A sample is a tree of cells. The tree depth depends on two things:
            - the difference between start and stop times
            - the average cell cycle time (set in DivisionParams)
        The number of evaluated values for an exaclty sampled process (e.g.
        the Ornstein Uhlenbeck process) is inversely proportional to the
        period at which data is recorded.

        The total number of operations depends:
            - exponentially on the difference 'stop - start'
            - inversely proportionally to 'period'
            - proportionally on nbr_container
            - proportionally on nbr_colony_per_container
    """

    def __init__(
        self,
        nbr_container=1,
        nbr_colony_per_container=1,
        start=0.0,
        stop=100.0,
        period=5.0,
        seed=None,
    ):
        self.nbr_container = nbr_container
        self.nbr_colony_per_container = nbr_colony_per_container
        self.start = start
        self.stop = stop
        self.period = period
        self.seed = seed
        # metadata helper
        self.content = [
            ("nbr_container", nbr_container),
            ("nbr_colony_per_container", nbr_colony_per_container),
            ("start", start),
            ("stop", stop),
            ("period", period),
            ("seed", seed),
        ]
        return

    def __str__(self):
        msg = "Simulation parameters:\n" "----------------------\n"
        left_column_size = 25
        right_column_size = 10
        for key, val in self.content:
            msg += "{}".format(key).ljust(left_column_size)
            msg += ": " + "{}\n".format(val).rjust(right_column_size)
        return msg


class DivisionParams(object):
    """Parameterize when a given cell divides depending on its observables
    
    This settings is based on the linear results explained elsewhere.
    
    Parameters
    ----------
    div_lambda : float between 0 and 1
        the lambda parameter is the mixing parameter between timer (0)
        and sizer (1) models. The adder model is obtained with lambda_ 0.5.
    div_size_target : float
        the size target sets the cutoff size for division. With *div_lambda=1*, and
        no fluctuations, division is triggered when size reaches exactly the size target
    div_mode : str {'gamma', 'fixed'}
        'fixed': only the linear response is taken;
        'gamma': adds fluctuations over linear response
    div_sd_to_mean : float (default 0.1)
        sets the standard deviation relative to mean value (default 10%) when used
        with the *div_mode='gamma'* option
    use_growth_rate : str {'parameter', 'birth'}
        division control scales as cell size growth rate, that can be given as
        a parameter (fixed), or from the actual growth rate value **at cell birth**
    """

    def __init__(
        self,
        div_lambda=0,
        div_size_target=2.0,
        div_mode="gamma",
        div_sd_to_mean=0.1,
        use_growth_rate="parameter",
    ):
        self.div_lambda = div_lambda
        self.div_size_target = div_size_target
        self.div_mode = div_mode
        self.div_sd_to_mean = div_sd_to_mean
        self.use_growth_rate = use_growth_rate
        self.content = [
            ("div_lambda", div_lambda),
            ("div_size_target", div_size_target),
            ("div_mode", div_mode),
            ("div_sd_to_mean", div_sd_to_mean),
            ("use_growth_rate", use_growth_rate),
        ]
        return

    def rv(self, birth_size=1.0, growth_rate=1.0 / 60 * np.log(2.0)):
        """Returns a random variate for interdivision time given *birth_size* and *growth_rate*.

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
        tau = (
            (1.0 - self.div_lambda) * np.log(2.0)
            + self.div_lambda * np.log(self.div_size_target / birth_size)
        ) / growth_rate
        if self.div_mode == "gamma":
            tau = _gamma(mean=tau, std=self.div_sd_to_mean * tau)
        return tau

    def __str__(self):
        msg = "Division parameters:\n" "--------------------\n"
        left_column_size = 5
        right_column_size = 10
        for key, val in self.content:
            msg = "{}".format(key).ljust(left_column_size)
            msg += ": " + "{}\n".format(val).rjust(right_column_size)
        return msg


def _gamma(mean=60, std=6.0):
    theta = std ** 2 / mean
    k = (std / theta) ** 2 + 1.0
    val = np.random.gamma(k, scale=theta)
    return val


class SampleInitialSize(object):
    """Initialize cell size given a target value for division and fluctuations.

    Parameters
    ----------
    birth_size_mean : float
    birth_size_mode : str {'fixed', 'lognormal'}
        initial distribution of birth size; fixed corresponds to a single value
        which is half the cutoff, lognormal returns a lognormal variate with
        parameter sigma
    birth_size_sd_to_mean : float
        relative fluctuations
        
    """

    def __init__(
        self, birth_size_mean=1.0, birth_size_mode="fixed", birth_size_sd_to_mean=0.1
    ):
        self.birth_size_mean = birth_size_mean
        self.birth_size_mode = birth_size_mode
        self.birth_size_sd_to_mean = birth_size_sd_to_mean
        self.content = [
            ("birth_size_mean", birth_size_mean),
            ("birth_size_mode", birth_size_mode),
            ("birth_size_sd_to_mean", birth_size_sd_to_mean),
        ]
        return

    def rv(self):
        if self.birth_size_mode == "fixed":
            return self.birth_size_mean
        elif self.birth_size_mode == "lognormal":
            # getting normal parameters from average and relative fluct
            av = self.birth_size_mean
            sigma = np.sqrt(np.log(1.0 + self.birth_size_sd_to_mean ** 2))
            mu = np.log(av) - 0.5 * sigma ** 2
            lognorm = np.random.lognormal(mean=mu, sigma=sigma)
            return lognorm


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

    def __init__(self, identifier=None, parent=None, birth_time=None, lifetime=None):

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
