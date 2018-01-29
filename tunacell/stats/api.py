#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module sets up api functions for dynamical correlation analysis.
"""
from __future__ import print_function

import numpy as np
import collections
import logging
import time

from tunacell.stats.utils import (iter_timeseries_,
                              iter_timeseries_2,
                              Region, Regions, UndefinedRegion,
                              CompuParams)
from tunacell.stats.single import (Univariate, StationaryUnivariate,
                               UnivariateIOError, StationaryUnivariateIOError)
from tunacell.stats.two import Bivariate, StationaryBivariate, StationaryBivariateIOError
from tunacell.stats.compute import (set_dynamics,
                                set_stationary_autocorrelation,
                                set_crosscorrelation,
                                set_stationary_crosscorrelation,
                                NoValidTimes)


logger = logging.getLogger(__name__)

MIN_INTERDIVISION_TIME = 5.  # World record is set by Vibrio natriegens
# See R.G.Eagon, J. Bact. vol 83, pp 736-737 (1962)


# SINGLE DYNAMIC ONBSERVABLE

def compute_univariate(exp, obs, region='ALL', cset=[], times=None,
                       size=None):
    """Computes one-point and two-point functions of statistical analysis.

    This functions handles conditions and time-window binning:

    * all conditions provided in cset are applied independently, in addition
      to the computation with unconditioned data (labelled 'master')
    * A time-binning window is provided with a given offset and a period.
      Explicitely a given time value t found in data will be rounded up to
      closest offset_t + k * delta_t, where k is an integer.

    Parameters
    ----------
    exp : :class:`Experiment` instance
    obs : :class:`Observable` instance
    region : :class:`Region` instance or str (default 'ALL')
        in case of str, must be the name of a registered region
    cset : list of :class:`FilterSet` instances
    times : 1d ndarray, or str (default None)
        array of times at which process is evaluated. Default is to use the
        'ALL' region with the period taken from experiment metadata. User can
        opt for a specific time array, or for the label of a region as a string
    size : int (default None)
        limit the iterator to size Lineage instances (used for testing)

    Returns
    -------
    Univariate instance
    """
    reg = _convert_region(region, exp)
    if isinstance(times, np.ndarray):
        eval_times = times
    elif isinstance(times, collections.Iterable):
        eval_times = np.array(times)
    else:
        eval_times = _default_eval_times(exp, obs, reg)
    # initialize Univariate and each of its item
    univ = Univariate(exp, obs, eval_times, reg, cset)  # empty
    # Set iterator over TimeSeries
    timeseries = iter_timeseries_(exp, obs, cset, size=size)
    # call the master function performing computation
    start_time = time.time()
    set_dynamics(timeseries, univ, eval_times)
    delta_time = time.time() - start_time
    msg = ('Univariate statistics for obs "{}"'.format(obs.name) + ''
           ', region "{}"'.format(reg.name) + ''
           'computed in {:.3f}'.format(delta_time))
    logger.info(msg)
    return univ


def _convert_region(region, exp):
    """Convert region when a string is used"""
    if isinstance(region, str):
        regs = Regions(exp)
        return regs.get(region)
    elif isinstance(region, Region):
        return region
    else:
        raise ValueError(region)


def _default_eval_times(exp, obs, region):
    """Returns default evalutation times.

    Parameters
    ----------
    exp : :class:`Experiment` instance
    obs : :class:`Observable` or :class:`FunctionalObservable` instance
    region : :class:`Region` instance
    """
    if obs.timing != 'g':
        period = exp.period
        tmin = region.tmin
        tmax = region.tmax
    else:
        period = 1
        n_max = int(np.ceil((region.tmax - region.tmin)/MIN_INTERDIVISION_TIME))
        tmin = - n_max
        tmax = n_max
    eval_times = np.arange(tmin, tmax + period, period)
    return eval_times


def load_univariate(exp, obs, region='ALL', cset=[]):
    """Initialize an empty Univariate instance.

    Such a Univariate instance is bound to an experiment (through exp),
    an observable, and a set of conditions.

    Parameters
    ----------
    exp : :class:`Experiment` instance
    obs : :class:`Observable` instance
    region : :class:`Region` instance or str (default 'ALL')
        in case of str, must be the name of a registered region
    cset : sequence of :class:`FilterSet` instances

    Returns
    -------
    :class:`Univariate` instance
        initialized, nothing computed yet

    Raises
    ------
    UnivariateIOError
        when importing fails (no data corresponds to input params)
    """
    logger.debug('Converting region and setting evaluation times')
    reg = _convert_region(region, exp)
    # use default eval_times to respect __init__, will be updated upon reading
    eval_times = _default_eval_times(exp, obs, reg)
    logger.debug('Instantiating univariate object for "{}"'.format(obs.name))
    univ = Univariate(exp, obs, eval_times, reg, cset)
    logger.debug('Trying to import results from files...')
    univ.import_from_text()
    logger.info('Import univariate statistics for "{}" successful'.format(obs.name))
    return univ


# SINGLE STATIONARY OBSERVABLE
    
class ParamsError(ValueError):
    pass


def _check_params(region, options):
    if not isinstance(options, CompuParams):
        raise ParamsError('options must be a CompuParams instance')
    if ((not hasattr(region, 'tmin')) or
         (not hasattr(region, 'tmax')) or
         (not hasattr(region, 'name'))):
        raise ParamsError('region must have name, tmin, tmax attributes')
    return


def load_stationary(univ, region, options):
    """Initialize a StationaryUnivariate instance from its dynamical one.

    Parameters
    -----------
    univ : :class:`Univariate` instance
    region : :class:`Region` instance
    options : :class:`CompuParams` instance

    Returns
    -------
    :class:`StationaryInstance` instance
        set up with empty arrays
    """
    _check_params(region, options)
    logger.debug('Instantiating stationary univariate object for "{}"'.format(univ.obs.name))
    stat = StationaryUnivariate(univ, region, options)
    logger.debug('Trying to import results from files...')
    stat.import_from_text()
    logger.info('Import stationary univariate statistics for "{}" successful'.format(univ.obs.name))
    _update_univariate_from_stationary(univ, stat)
    return stat


def _update_univariate_from_stationary(univ, stat):
    for lab in stat._condition_labels:
        cdt_stat = stat[lab]
        univ[lab].stationary = cdt_stat
    return


def compute_stationary(univ, region, options, size=None):
    """Computes stationary autocorrelation. API level.

    Parameters
    ----------
    univ : :class:`Univariate` instance
        the stationary autocorr is based on this object
    region : :class:`Region` instance
    options : :class:`CompuParams` instance
        set the 'adjust_mean' and 'disjoint' options
    size : int (default None)
        limit number of parsed Lineages

    """
    _check_params(region, options)
    # initialize StationaryUnivariate
    stationary = StationaryUnivariate(univ, region, options)
    # Set iterator over TimeSeries
    timeseries = iter_timeseries_(univ.exp, univ.obs, univ.cset, size=size)
    # call the function performing computation and updating stationary
    try:
        start_time = time.time()
        set_stationary_autocorrelation(timeseries, univ, stationary,
                                       tmin=region.tmin, tmax=region.tmax,
                                       adjust_mean=options.adjust_mean,
                                       disjoint=options.disjoint)
        delta_time = time.time() - start_time
        msg = ('Steady statistics for "{}"'.format(univ.obs.name) + ' '
               'over region "{}"'.format(region.name) + ' '
               'computed in {:.3f}'.format(delta_time))
        logger.info(msg)
        _update_univariate_from_stationary(univ, stationary)
    except NoValidTimes as error:
        print('No valid time points for {} in {}'.format(univ.obs.name, region))
        raise error
    return stationary


# CROSS CORRELATION
    
def _update_univariate_from_bivariate(univs, two):
    for cdt_lab in two._condition_labels:
        cdt_two = two[cdt_lab]
        for k in range(2):
            univs[k][cdt_lab].two[str(univs[k-1].obs)] = cdt_two
    return


def load_bivariate(row_univariate, col_univariate):
    """Initialize a StationaryBivariate instance from its dynamical one.

    Parameters
    -----------
    row_univariate : :class:`Univariate` instance
    col_univariate : :class:`Univariate` instance

    Returns
    -------
    :class:`Bivariate` instance
        set up with empty arrays
    """
    two = Bivariate(row_univariate, col_univariate)
    two.import_from_text()
    univs = row_univariate, col_univariate
    _update_univariate_from_bivariate(univs, two)
    return two


def compute_bivariate(row_univariate, col_univariate, size=None):
    """Computes cross-correlation between observables defiend in univs.

    This functions handles conditions and time-window binning:

    * all conditions provided in cset are applied independently, in addition
      to the computation with unconditioned data (labelled 'master')
    * A time-binning window is provided with a given offset and a period.
      Explicitely a given time value t found in data will be rounded up to
      closest offset_t + k * delta_t, where k is an integer.

    Parameters
    ----------
    univs : couple of Univariate instances
    size : int (default None)
        limit the iterator to size Lineage instances (used for testing)

    Returns
    -------
    TwoObservable instance
    """
    univs = row_univariate, col_univariate
    s1, s2 = univs
    obs1 = s1.obs
    obs2 = s2.obs
    # check that the two observables are either two dynamics mode, either zero
    if obs1.mode != obs2.mode:
        if obs1.mode == 'dynamics' or obs2.mode == 'dynamics':
            msg = ('Cannot mix time-lapse and cell-cycle observables:\n'
                   'Here obs1.mode = {} and obs2.mode = {}'.format(obs1.mode,
                                                                   obs2.mode))
            raise TypeError(msg)
        elif obs1.timing == 'g':
            if obs2.timing != 'g':
                msg = ('If one observable is evaluated in generation time, '
                       'the other one must be as well in generation time.')
                raise TypeError(msg)
        elif obs2.timing == 'g':
            if obs1.timing != 'g':
                msg = ('If one observable is evaluated in generation time, '
                       'the other one must be as well in generation time.')
                raise TypeError(msg)
    # initialize Univariate and each of its item
    two = Bivariate(row_univariate, col_univariate)  # empty
    exp = two.exp
    cset = two.cset
    start_time = time.time()
    timeseries = iter_timeseries_2(exp, obs1, obs2, cset, size=size)
    # call the master function performing computation
    set_crosscorrelation(timeseries, row_univariate, col_univariate, two)
    delta_time = time.time() - start_time
    msg = ('Bivariate statistics for "{}--{}"'.format(obs1.name, obs2.name) + ' '
           'computed in {:.3f}'.format(delta_time))
    logger.info(msg)
    # update conditioned univ cross-correlation
    _update_univariate_from_bivariate(univs, two)
    return two


def _update_univariate_from_stationary_bivariate(univs, stwo):
    for cdt_lab in stwo._condition_labels:
        cdt_two = stwo[cdt_lab]
        for k in range(2):
            univs[k][cdt_lab].two_stationary[str(univs[k-1].obs)] = cdt_two
    return


def load_stationary_bivariate(row_univariate, col_univariate,
                              region, options):
    """Initialize a StationaryBivariate instance from its dynamical one.

    Parameters
    -----------
    row_univariate : :class:`Univariate` instance
    col_univariate : :class:`Univariate` instance
    region : :class:`Region` instance
    options : :class:`CompuParams` instance

    Returns
    -------
    :class:`StationaryBivariate` instance
        set up with empty arrays
    """
    stwo = StationaryBivariate(row_univariate, col_univariate,
                               region, options)
    stwo.import_from_text()
    univs = row_univariate, col_univariate
    # stationary univariate need to be computed for full inspection
    for univ in univs:
        try:
            suniv = load_stationary(univ, region, options)
        except StationaryUnivariateIOError:
            suniv = compute_stationary(univ, region, options)
            suniv.export_text()
        _update_univariate_from_stationary(univ, suniv)
    _update_univariate_from_stationary_bivariate(univs, stwo)
    return stwo


def compute_stationary_bivariate(row_univariate, col_univariate,
                                 region, options, size=None):
    """Computes stationary cross-correlation function from couple of univs

    Need to compute stationary univariates as well.
    """
    s1, s2 = row_univariate, col_univariate
    obs1 = s1.obs
    obs2 = s2.obs
    univs = row_univariate, col_univariate
    # stationary univariate need to be computed for full inspection
    for univ in univs:
        try:
            suniv = load_stationary(univ, region, options)
        except StationaryUnivariateIOError:
            suniv = compute_stationary(univ, region, options)
            suniv.export_text()
        _update_univariate_from_stationary(univ, suniv)
    sbivar = StationaryBivariate(row_univariate, col_univariate,
                                 region, options)
    exp = sbivar.exp
    cset = sbivar.cset
    start_time = time.time()
    timeseries = iter_timeseries_2(exp, obs1, obs2, cset, size=size)
    set_stationary_crosscorrelation(timeseries, row_univariate, col_univariate,
                                    sbivar,
                                    tmin=region.tmin, tmax=region.tmax,
                                    adjust_mean=options.adjust_mean,
                                    disjoint=options.disjoint)
    delta_time = time.time() - start_time
    msg = ('Steady bivariate statistics for "{}--{}"'.format(obs1.name, obs2.name) + ' '
           'over region "{}"'.format(region.name) + ' '
           'computed in {:.3f}'.format(delta_time))
    logger.info(msg)
    # update conditioned univ stationary cross-correlation
    _update_univariate_from_stationary_bivariate(univs, sbivar)
    return sbivar
