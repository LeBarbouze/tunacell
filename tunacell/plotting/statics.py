#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module is deprecated...
"""
from __future__ import print_function

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from scipy import stats
from scipy.stats import spearmanr
from matplotlib.patches import Polygon

from tunacell.base.datatools import gaussian_smooth

def scatter_plot(bivariate, xsize=2.3, ysize=2.1,
                 xrange=(None, None),
                 use_xname=None,
                 yrange=(None, None),
                 use_yname=None,
                 groupby=None,
                 which_keys='all',
                 color_index=0,
                 bins='auto', xunits='', yunits=''):
    """Draw a scatter plot and empirical distributions
    
    Parameters
    ----------
    bivariate : StationaryBivariate instance
    xsize : float
        horizontal size of plots (inches)
    ysize : float
        vertical size of plots(inches)
    xrange : tuple (float, float)
    use_xname : str (default None)
        x-observable name to use (default uses observable.name)
    yrange : tuple (float, float)
    use_yname : str (default None)
        y-observable name to use
    groupby : str
        pandas groupby technique to group values according to a key; available
        keys are condition representations ('FilterSet(...)') or 'g' when
        generation timing is used
    which_keys : list, or str (default 'all')
        when groupby is used, which key value to plot ('True' for example in
        case of condition, or 1 for plotting only generation 1)
    color_index : int (default 0)
        set first color to be used in matplotlib color cycle
    bins : int, or str
        see numpy.histogram
    xunits : str
    yunits : str
    
    Returns
    -------
    fig : Figure instance
    ax0 : Axes
        scatter plot axes
    ax1 : Axes
        x-distribution axes
    ax2 : Axes
        y-distribution axes
    """
    fig = plt.figure(figsize=(xsize, ysize))
    gs = gridspec.GridSpec(3, 3)
    gs.update(wspace=0.0, hspace=0.0)
    ax1 = fig.add_subplot(gs[0, :2])  # x distribution
    ax0 = fig.add_subplot(gs[1:, :2])  # scatter plot
    ax2 = fig.add_subplot(gs[1:, 2])  # y distribution
    ax1.tick_params(axis='x', labelbottom='off', top='on', labeltop='on')
    ax1.tick_params(axis='y', labelleft='off', left='off')
    ax2.tick_params(axis='y', labelleft='off', labelright='on', right='on')
    ax2.tick_params(axis='x', labelbottom='off', bottom='off')
    ax0.tick_params(axis='x', top='on', direction='in')

    u1, u2 = bivariate.univariates
    o1, o2 = u1.obs, u2.obs
    
    handles = []
    xs = []
    ys = []

    df = bivariate.dataframe
    if groupby is None:
        x_ = df[o1.name].values
        y_ = df[o2.name].values
    
        data, = ax0.plot(x_, y_, ls='None', marker='o', alpha=.3,
                         color='C{}'.format(color_index),
                         label='data ({} samples)'.format(len(x_)))
        handles.append(data)
        xs.append(x_)
        ys.append(y_)
        color_index += 1
    else:
        gs = df.groupby(groupby)
        if which_keys == 'all':
            keys = gs.groups.keys()
        else:
            keys = [key for key in gs.groups.keys() if key in which_keys]
        for key in keys:
            index = gs.groups[key]
            x_ = df[o1.name].loc[index].values
            y_ = df[o2.name].loc[index].values
            
            data, = ax0.plot(x_, y_, ls='None', marker='o', alpha=.3,
                             color='C{}'.format(color_index),
                             label='{}={} data ({} samples)'.format(groupby, key, len(x_)))
            handles.append(data)
            xs.append(x_)
            ys.append(y_)
            color_index += 1

    if xrange[0] is not None:
        for ax in [ax0, ax1]:
            ax.set_xlim(left=xrange[0])
    if xrange[1] is not None:
        for ax in [ax0, ax1]:
            ax.set_xlim(right=xrange[1])
    if yrange[0] is not None:
        for ax in [ax0, ax2]:
            ax.set_ylim(bottom=yrange[0])
    if yrange[1] is not None:
        for ax in [ax0, ax2]:
            ax.set_ylim(top=yrange[1])
    
    for x, y, h in zip(xs, ys, handles):
        color = h.get_color()
        # x distribution
        left, right = ax0.get_xlim()
        xh, xbins = np.histogram(x, bins=bins, range=(left, right), density=True)
        xcoords = (xbins[:len(xbins)-1] + xbins[1:])/2.
        ax1.plot(xcoords, xh, color=color)
    
        # y distribution
        bottom, top = ax0.get_ylim()
        yh, ybins = np.histogram(y, bins=bins, range=(bottom, top), density=True)
        ycoords = (ybins[:len(ybins)-1] + ybins[1:])/2.
        ax2.plot(yh, ycoords, color=color)

    # labels
    xlabel = r'{}'.format(o1.latexify(show_variable=False, use_name=use_xname))
    if xunits:
        xlabel += ' ({})'.format(xunits)
    ylabel = r'{}'.format(o2.latexify(show_variable=False, use_name=use_yname))
    if yunits:
        ylabel += ' ({})'.format(yunits)
    ax0.set_xlabel(xlabel, size='medium')
    ax0.set_ylabel(ylabel, size='medium')

    return fig, ax0, ax1, ax2, handles
