#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module is deprecated...
"""

import numpy as np

import matplotlib.pyplot as plt

from scipy import stats
from scipy.stats import spearmanr
from matplotlib.patches import Polygon

from tuna.plotting.defs import params, colors
from tuna.base.datatools import gaussian_smooth


def scatter(dat, xkey=None, ykey=None, delta=1000, xlim=[], ylim=[],
            thresh=[0, 10, 50, 90, 100],
            sig=0.06,
            jupyter=False,
            savefig=None):
    if xkey is None or ykey is None:
        print "No keys"
        return
    elif (xkey not in dat.dtype.fields.keys() or
          ykey not in dat.dtype.fields.keys()):
        print "Unvalid keys"
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
    print "there's {0} valid samples".format(len(so))

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
    print percs
    slices = [(start, stop) for start, stop in zip(percs[:-1], percs[1:])]

    xx = np.linspace(ybottom, ytop, 100)

    for i, (start, stop) in enumerate(slices):
        color = colors[i % len(colors)]
        vals = ys[np.logical_and(np.greater(xs, start),
                                 np.less_equal(xs, stop))]
        vals = np.sort(vals)
        lab = '{0}: {1} in ({2:.3f}, {3:.3f}), \
              '.format(i+1, xkey, slices[i][0], slices[i][1])
        print lab
        print 'm: {0:.3f}, std:{1:.3f}, skew: {2:.3f}\
              '.format(np.mean(vals), np.std(vals), stats.skew(vals))
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
