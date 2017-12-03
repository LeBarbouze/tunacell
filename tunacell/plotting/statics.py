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


def scatter(dat, xkey=None, ykey=None, delta=1000, xlim=[], ylim=[],
            thresh=[0, 10, 50, 90, 100],
            sig=0.06,
            jupyter=False,
            savefig=None):
    if xkey is None or ykey is None:
        print("No keys")
        return
    elif (xkey not in dat.dtype.fields.keys() or
          ykey not in dat.dtype.fields.keys()):
        print("Unvalid keys")
        return

    x, y = map(np.array, zip(*dat[[xkey, ykey]]))
    # convert rate from /min to /hour
    if 'rate' in xkey:
        x = x * 60.
    if 'rate' in ykey:
        y = y * 60.

    couples = zip(x, y)

    # Pearson correlation coefficient and p-value
    pr = spearmanr(*zip(*couples))
    so = sorted(couples, key=lambda item: item[0])
    so = [co for co in so if not np.isnan(co[1])]
    print("there's {0} valid samples".format(len(so)))

    # delta = 1000
    cuts = np.arange(0, len(so), delta)
    sls = [slice(cuts[i], cuts[i+1]) for i in range(0, len(cuts)-1)] + \
          [slice(cuts[-1], len(so))]
    sxx, syy = zip(*so)
    binned = [(np.nanmean(sxx[sl]), np.nanmean(syy[sl]), np.nanstd(syy[sl]))
              for sl in sls]

    if xlim:
        xleft, xright = xlim
    else:
        xleft, xright = np.amin(sxx), np.amax(sxx)
    if ylim:
        ybottom, ytop = ylim
    else:
        ybottom, ytop = np.amin(syy), np.amax(syy)

    fig, axs = plt.subplots(2, 2, figsize=(16, 12))

    ax = axs[0, 0]

    ax.scatter(x, y, facecolor='b', alpha=0.4, edgecolor='none')

    ax.text(0.6, 0.9, 'Spearman r: {0:.2f}'.format(pr[0]),
            transform=ax.transAxes, fontsize=16)

    if xkey in params.keys():
        if 'unit' in params[xkey].keys():
            xlabel = r''+xkey+' ({0})'.format(params[xkey]['unit'])
        else:
            xlabel = xkey
    else:
        xlabel = xkey
    if ykey in params.keys():
        if 'unit' in params[ykey].keys():
            ylabel = r''+ykey+' ({0})'.format(params[ykey]['unit'])
        else:
            ylabel = ykey
    else:
        ylabel = ykey

    ax.set_xlabel(xlabel, size=16)
    ax.set_ylabel(ylabel, size=16)

    ax.set_xlim(left=xleft, right=xright)
    ax.set_ylim(bottom=ybottom, top=ytop)

    bins, vals, stds = map(np.array, zip(*binned))
    plus = vals + 0.5*stds
    minus = vals - 0.5*stds

    points = zip(bins[::-1], minus[::-1]) + zip(bins, plus)
    cpoints = [pt for pt in points if not np.isnan(pt[1])]

    pol = Polygon(cpoints, closed=True, color='r', alpha=0.4)
    ax.add_patch(pol)
    ax.plot(bins, vals, marker='o',
            markersize=8.,
            markeredgecolor='r',
            markerfacecolor='None',
            ls='-', lw=2, color='r')

    vals = y
    h, bins = np.histogram(vals, bins=20, range=(ybottom, ytop), normed=True)
    coords = (bins[:len(bins)-1] + bins[1:])/2.
    axs[0, 1].plot(coords, h,  marker='o',
                   ms=8.,
                   ls='none',
                   markeredgecolor='b',
                   markeredgewidth=1.,
                   markerfacecolor='b',
                   alpha=0.6)
    xx = np.linspace(ybottom, ytop, 100)
    x, yy = zip(*gaussian_smooth(coords, h, sigma=sig, x=xx))

    axs[0, 1].plot(x, yy, '--', lw=2, color='b')
    axs[0, 1].set_xlim(left=ybottom, right=ytop)
    axs[0, 1].set_xlabel(ylabel, size=16)
    axs[0, 1].set_ylabel('Probability distribution', size=16)

    xs, ys = map(np.array, zip(*so))

    thresh = [0, 10, 50, 90, 100]
    percs = map(lambda th: np.percentile(xs, th), thresh)
    print(percs)
    slices = [(start, stop) for start, stop in zip(percs[:-1], percs[1:])]

    xx = np.linspace(ybottom, ytop, 100)

    for i, (start, stop) in enumerate(slices):
        color = colors[i % len(colors)]
        vals = ys[np.logical_and(np.greater(xs, start),
                                 np.less_equal(xs, stop))]
        vals = np.sort(vals)
        lab = '{0}: {1} in ({2:.3f}, {3:.3f}), \
              '.format(i+1, xkey, slices[i][0], slices[i][1])
        print(lab)
        print('m: {0:.3f}, std:{1:.3f}, skew: {2:.3f}\
              '.format(np.mean(vals), np.std(vals), stats.skew(vals)))
        axs[1, 0].plot(vals, np.array(range(len(vals)))/float(len(vals)),
                       ls='-', color=color, lw=2,
                       label=lab)

        h, bins = np.histogram(vals, bins=20, range=(ybottom, ytop),
                               normed=True)
        coords = (bins[:len(bins)-1] + bins[1:])/2.
        axs[1, 1].plot(coords, h,  marker='o',
                       ms=8.,
                       ls='none',
                       markeredgecolor=color,
                       markeredgewidth=1.,
                       markerfacecolor=color,
                       alpha=0.6,
                       label='{0}'.format(i+1))

        x, yy = zip(*gaussian_smooth(coords, h, sigma=sig, x=xx))

        axs[1, 1].plot(x, yy, '--', lw=2, color=color,
                       label='{0}'.format(i+1))

    axs[1, 0].legend(loc=4)
    axs[1, 0].set_ylabel('Cumulative distribution', size=16)
    axs[1, 1].set_ylabel('Distribution (Gaussian smoothed)', size=16)
    for ax in axs[1, :]:
        ax.set_xlim(left=ybottom, right=ytop)
        ax.set_xlabel(ylabel, size=16)

    if not jupyter:
        plt.show()

    if savefig is not None:
        fig.savefig(savefig)

    return
