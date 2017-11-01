#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Helper functions for plotting modules
"""

import numpy as np
from matplotlib import ticker


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



def _set_timelabel(obs):
    """For a given observable, returns the timelabel to be used in plots"""
    if obs.mode != 'dynamics' and obs.timing == 'g':
        timelabel = 'Generations'
        if obs.tref is not None:
            timelabel += ' (since tref {})'.format(obs.tref)
    else:
        timelabel = 'Time (minutes)'
    return timelabel