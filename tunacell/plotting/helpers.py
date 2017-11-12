#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Helper functions for plotting modules
"""
import collections

import numpy as np
from matplotlib import ticker

from tunacell.base import Colony, Lineage

MIN_SEP = 20
MAX_XTICKS = 6


def _set_axis_limits(ax, values, which='x', pad=.1, force_range=(None, None)):
    """Computes axis limits with padding given values.

    Parameters
    ----------
    ax : Axes instance
        axes tp set limits
    values : iterable/list of ints/floats
        list of values to find bounds
    which : str {'x', 'y'}
        which axis to set limits
    pad : float (default .1)
        fractional padding to add at both boundaries
    force_range: couple of floats (default (None, None))
        when not None, will force this boun

    Raises
    ------
    ValueError
        when which is not 'x' or 'y'
    """
    if len(values) == 0:
        lower, upper = -1, 1  # default
    else:
        lower, upper = np.nanmin(values), np.nanmax(values)
    if force_range[0] is not None:
        lower = force_range[0]
    if force_range[1] is not None:
        upper = force_range[1]
    if lower == upper:
        vrange = .1 * upper  # 10% of unique value
    else:
        vrange = upper - lower
    low = lower - pad * vrange
    up = upper + pad * vrange
    if which == 'x':
        ax.set_xlim(left=low, right=up)
    elif which == 'y':
        ax.set_ylim(bottom=low, top=up)
    else:
        raise ValueError("'which' argument can be either 'x' or 'y'")
    return lower, upper


def _set_time_axis_ticks(ax, obs, bounds=(None, None)):
    """Set ticker options for time axis
    
    Parameters
    ----------
    ax : Axes instance
    obs : Observable instance
    bounds : couple of floats
        left, right values
    
    Returns
    -------
    ticker locator
    """
    left, right = bounds
    found = ax.get_xlim()
    if left is None:
        left = found[0]
    if right is None:
        right = found[1]
    locator = ticker.AutoLocator()
    if obs.timing == 'g':
        locator = ticker.MaxNLocator(bins='auto', integer=True)
    else:
        # set lowest tick sep to MIN_SEP
        n_seps = 1000  # large number to initialize
        multiple = 0
        while n_seps > MAX_XTICKS and multiple < 300:  # later correspond to 100 hours sep between ticks... unlikely
            multiple += 1
            sep = multiple * MIN_SEP
            n_seps = int(np.ceil((right - left) / sep))
        if 2 <= n_seps <= MAX_XTICKS:
            # for range > 1 hour, get only a multiple of hours
            if sep > 60:
                hour_multiple = int(np.ceil(sep / 60))
                sep = hour_multiple * 60
            locator = ticker.MultipleLocator(base=sep)
    return locator


def _set_timelabel(obs, use_tref=True):
    """For a given observable, returns the timelabel to be used in plots

    Parameters
    ----------
    obs : Observable instance
    use_tref : bool
        whether to use tref (set False for stationary plotting)

    Returns
    -------
    timelabel : str
    """
    if obs.mode != 'dynamics' and obs.timing == 'g':
        timelabel = 'Generations'
        if obs.tref is not None and use_tref:
            timelabel += ' (since tref {})'.format(obs.tref)
    else:
        timelabel = 'Time (minutes)'
    return timelabel


def _unroll_samples(arg):
    """Iterator over Lineage instances"""
    if isinstance(arg, Lineage):
        yield arg
    elif isinstance(arg, Colony):
        for lin in arg.iter_lineages():
            yield lin
    elif isinstance(arg, collections.Iterable):
        for item in arg:
            for elem in _unroll_samples(item):
                yield elem