#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module defines plotting functions for the statistics of the dynamics.
"""
from __future__ import print_function

import os
import warnings
import numpy as np
import string
import collections
import logging

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, ScalarFormatter
import matplotlib.transforms as transforms
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker

from tunacell.filters.main import FilterSet
from tunacell.stats.single import Univariate, StationaryUnivariate
from tunacell.stats.two import StationaryBivariate
from tunacell.io import text


class UnivariatePlot(object):
    """Class to make onepoint and twopoint plotting

    .. note:: Deprecated in tunacell 0.0.7

    Parameters
    ----------
    univariate : :class:`Univariate` instance
    """

    def __init__(self, univariate):
        self.univariate = univariate
        self.fig1 = None  # one-point plot
        self.fig2 = None  # two-point plot
        return

    def make_onepoint(self, show_cdts='all', left=None, right=None,
                      show_ci=True,
                      mean_ref=None, mean_lims=(None, None),
                      var_ref=None, var_lims=(None, None)):
        fig = plot_onepoint(self.univariate, show_cdts=show_cdts,
                            left=left, right=right,
                            mean_ref=mean_ref, mean_lims=mean_lims,
                            show_ci=show_ci,
                            var_ref=var_ref, var_lims=var_lims)
        self.fig1 = fig
        return

    def make_twopoints(self, condition_label=None, trefs=[], ntrefs=3,
                       trange=(-100., 100.)):
        fig = plot_twopoints(self.univariate, condition_label=condition_label,
                             trefs=trefs, ntrefs=ntrefs, trange=trange)
        self.fig2 = fig
        return

    def save(self, user_path=None, force_export_text=False,
             label='plot', extension='.pdf'):
        # one point plot is saved under master
        univ = self.univariate
        try:
            obs_path = univ._get_obs_path(user_root=user_path,
                                            write=False)
        except text.MissingFolderError:
            # it means data has not been written yet
            # export data and then get
            univ.export_text()
            obs_path = univ._get_obs_path(user_root=user_path,
                                            write=False)
        # export text jointly
        if force_export_text:
            univ.export_text(analysis_folder_path=user_path)
        if self.fig1 is not None:
            name = ''
            if label:
                name += label + '_'
            name += 'onepoint' + extension
            figname = os.path.join(obs_path, name)
            self.fig1.savefig(figname, bbox_inches='tight')
        else:
            print('Use make_onepoint(), then save')
        if self.fig2 is not None:
            name = ''
            if label:
                name += label + '_'
            name += 'twopoints' + extension
            figname = os.path.join(obs_path, name)
            self.fig2.savefig(figname, bbox_inches='tight')
        else:
            print('Use make_twopoints(), then save')
        return


def _append_cdt(univariate, this_cdt, cdt_list):
    """Append condition associated to this_cdt in univariate object to cdt_list

    Parameters
    ----------
    univariate : :class:`Univariate` instance
    this_cdt : str or :class:`FilterSet` instance
        either the condition instance or its string representation
    cdt_list : list of conditions
        list of conditions to append condition to
    """
    found = False
    if isinstance(this_cdt, str):
        # find which
        for cdt in univariate.cset:
            if repr(cdt) == this_cdt:
                found = True
                break
    elif isinstance(this_cdt, FilterSet):
        for cdt in univariate.cset:
            if repr(cdt) == repr(this_cdt):
                found = True
                break
    if found:
        cdt_list.append(cdt)
    return


def plot_onepoint(univariate, show_cdts='all', show_ci=False,
                  mean_ref=None, var_ref=None,
                  axe_xsize=6., axe_ysize=2.,
                  time_range=(None, None),
                  time_fractional_pad=.1,
                  counts_range=(None, None),
                  counts_fractional_pad=.1,
                  average_range=(None, None),  # auto
                  average_fractional_pad=.1,
                  variance_range=(None, None),
                  variance_fractional_pad=.1,
                  show_legend=True,
                  show_cdt_details_in_legend=False,
                  save=False, user_path=None, ext='.png'):
    """Plot one point statistics.

    Parameters
    ----------
    univariate : Univariate instance
    show_cdts : str (default 'all')
        must be either 'all', or 'master', or the repr of a condition, or a
        list thereof
    mean_ref : float
        reference mean value: what user expect to see as sample average to
        compare with data
    var_ref : float
        reference variance value: what user expect to see as sample variance to
        compare with data
    axe_xsize : float (default 6)
        size of the x-axis (inches)
    axe_ysize : float (default 2.)
        size if a single ax y-axis (inches)
    time_range : couple of floats (default (None, None))
        specifies (left, right) bounds
    time_fractional_pad : float (default .1)
        fraction of x-range to add as padding
    counts_range : couple of floats (default (None, None))
        specifies range for the Counts y-axis
    counts_fractional_pad : float (default .2)
        fractional amount of y-range to add as padding
    average_range : couple of floats (default (None, None))
        sepcifies range for the Average y-axis
    average_fractional_pad : couple of floats (default .2)
        fractional amounts of range to padding
    variance_range : couple of floats (default (None, None))
        sepcifies range for the Variance y-axis
    average_fractional_pad : couple of floats (default .2)
        fractional amounts of range to padding
    show_legend : bool {True, False}
        print out legend
    show_cdt_details_in_legend : bool {False, True}
        show details about filters
    save : bool {False, True}
        whether to save plot
    user_path : str (default None)
        user defined path where to save figure; default is canonical path
        (encouraged)
    ext : str {'.png', '.pdf'}
        extension to be used when saving file
    """
    if not isinstance(univariate, Univariate):
        raise TypeError('Input is not {}'.format(Univariate))

    fig, axs = plt.subplots(3, 1, figsize=(axe_xsize, 3*axe_ysize))

    obs = univariate.obs

    # define time label
    if obs.mode != 'dynamics' and obs.timing == 'g':
        timelabel = 'Generations'
        if univariate.obs.tref is not None:
            timelabel += ' (since tref {})'.format(univariate.obs.tref)
        # force ticks to be integers
        for ax in axs:
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    else:
        timelabel = 'Time (mins)'

    main_handles = []
    additional_handles = []

    all_times = []
    all_counts = []
    all_average = []
    all_variance = []

    # build condition list
    conditions = ['master', ]  # list of conditions to be plotted
    if show_cdts == 'all':
        conditions = ['master', ] + univariate.cset
    elif show_cdts == 'master':
        pass
    elif isinstance(show_cdts, collections.Iterable):
        for item in show_cdts:
            _append_cdt(univariate, item, conditions)
    else:
        _append_cdt(univariate, show_cdts, conditions)

    default_lw = mpl.rcParams['lines.linewidth']
    for index, cdt in enumerate(conditions):

        if cdt == 'master':
            c_repr = 'master'
            c_label = 'all samples'
            lw = default_lw + 1
            alpha = 1
            alpha_fill = .5
        else:
            c_repr = repr(cdt)
            if show_cdt_details_in_legend:
                c_label = str(cdt)
            else:
                c_label = cdt.label
            lw = default_lw
            alpha = .8
            alpha_fill = 0.3

        ok = np.where(univariate[c_repr].count_one > 0)

        times = univariate[c_repr].time[ok]
        all_times.extend(times)
        counts = univariate[c_repr].count_one[ok]
        all_counts.extend(counts)
        mean = univariate[c_repr].average[ok]
        all_average.extend(mean)
        var = univariate[c_repr].var[ok]
        all_variance.extend(var)
        std = univariate[c_repr].std[ok]
        se = 2.58 * std / np.sqrt(counts)
#        var = np.diagonal(univariate[c_repr].autocorr)

        line_counts, = axs[0].plot(times, counts, alpha=alpha, lw=lw,
                               label='{}'.format(c_label))
        main_handles.append(line_counts)
        color = line_counts.get_color()

        average, = axs[1].plot(times, mean, color=color, alpha=0.8, lw=lw, label=c_label)
        if show_ci:
            fill_std = axs[1].fill_between(times, mean-se, mean+se,
                                       facecolor=color, alpha=alpha_fill)
            # add onbly if empty (no need to repeat on conditions)
            additional_handles.append(fill_std)
            all_average.extend(mean-se)
            all_average.extend(mean+se)

        variance, = axs[2].plot(times, var, color=color, alpha=0.8, lw=lw, label=c_label)

    # adding reference lines
    if mean_ref is not None:
        mref = axs[1].axhline(mean_ref, ls='-.', color='C7', alpha=.7,
                               label='reference value')
        main_handles.append(mref)
        all_average.append(mean_ref)
    if var_ref is not None:
        vref = axs[2].axhline(var_ref, ls='-.', color='C7', alpha=.7,
                              label='reference value')
        # check last label if meanèref has been saved
        last_lab = main_handles[-1].get_label()
        if last_lab != vref.get_label():
            main_handles.append(vref)
        all_variance.append(var_ref)

    # xaxis limits
    left, right = np.nanmin(all_times), np.nanmax(all_times)
    if time_range[0] is not None:
        left = time_range[0]
    if time_range[1] is not None:
        right = time_range[1]
    hrange = right - left
    for ax in axs:
        ax.set_xlim(left=left - time_fractional_pad*hrange,
                    right=right + time_fractional_pad*hrange)

    # yaxis limits
    max_count = np.nanmax(all_counts)
    axs[0].set_ylim(bottom=0, top=max_count*(1. + counts_fractional_pad))  # counts

    # average
    bottom, top = np.nanmin(all_average), np.nanmax(all_average)
    if average_range[0] is not None:
        bottom = average_range[0]
    if average_range[1] is not None:
        top = average_range[1]
    vrange = top - bottom
    axs[1].set_ylim(bottom=bottom - average_fractional_pad*vrange,
                    top=top + average_fractional_pad*vrange)
    
    # variance
    bottom, top = np.nanmin(all_variance), np.nanmax(all_variance)
    if variance_range[0] is not None:
        bottom = variance_range[0]
    if average_range[1] is not None:
        top = variance_range[1]
    vrange = top - bottom
    axs[2].set_ylim(bottom=bottom - variance_fractional_pad*vrange,
                    top=top + variance_fractional_pad*vrange)

    # print vertical line at tref
    if univariate.obs.timing != 'g' and isinstance(univariate.obs.tref, float):
        for ax in axs:
            ax.axvline(univariate.obs.tref, color='C7', ls='--', alpha=.5)

    # ticks and labels
    for ax in axs:
        ax.xaxis.set_major_locator(ticker.MultipleLocator(_spacing_trange(hrange)))
    # counts
    axs[0].yaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
    for ax in axs[1:]:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=2))
    formatter = ScalarFormatter(useMathText=True, useOffset=False)
    formatter.set_powerlimits((-2, 4))
    for ax in axs:
        ax.yaxis.set_major_formatter(formatter)
        t = ax.yaxis.get_offset_text()
        plt.draw()
        msg = t.get_text()
        ax.text(0, .95, msg, ha='left', va='top', transform=ax.transAxes)
        t.set_visible(False)

    axs[0].tick_params(axis='x', direction='in', bottom='on', labelbottom='on', pad=-10)
    axs[1].tick_params(axis='x', direction='in', bottom='on', labelbottom='on', pad=-10)
    axs[2].set_xlabel(timelabel, x=.95, horizontalalignment='right',
                      fontsize='large')

    # hide intermediate x axis
    for ax in axs[:2]:
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='x', colors='C7')
    for ax in axs[1:]:
        ax.spines['top'].set_color('C7')

    axs[0].set_ylabel('Counts', fontsize='large')
    axs[1].set_ylabel('Average', fontsize='large')
    axs[2].set_ylabel('Variance', fontsize='large')
    # legend
    handles = main_handles[:]
    labels = [h.get_label() for h in handles]
    axs[-1].legend(handles=handles, labels=labels, loc='upper left',
                   bbox_to_anchor=(0, -.6/axe_ysize))
    # C.I.
    if additional_handles:
        ci = additional_handles[0]
        ci.set_color('C7')
        ci.set_label('.99 C.I.')
        labels = [item.get_label() for item in additional_handles]
        axs[1].legend(handles=[ci, ], labels=[ci.get_label(), ], loc=0)

    # add info above first ax
#    info = '(binning interval {})'.format(univariate.indexify.binsize)
#    axs[2].text(0., -.2, info, va='top', transform=axs[2].transAxes)

    # title
    latex_obs = obs.as_latex_string
    axs[0].text(0.5, 1+.2/axe_ysize,
                r'Statistics: {}'.format(latex_obs),
                size='x-large',
                horizontalalignment='center',
                verticalalignment='bottom',
                transform=axs[0].transAxes)
    fig.subplots_adjust(hspace=0)
    if save:
        univ = univariate
        try:
            obs_path = univ._get_obs_path(user_root=user_path, write=False)
        except text.MissingFolderError:
            # it means data has not been written yet
            # export data and then get
            univ.export_text(analysis_folder=user_path)
            obs_path = univ._get_obs_path(user_root=user_path, write=False)
        bname = 'plot_onepoint_' + univ.obs.name + '_' + univ.region.name + ext
        fname = os.path.join(obs_path, bname)
        fig.savefig(fname, bbox_to_inches='tight', pad_inches=0)
    return fig


def _spacing_trange(hrange):
    if hrange < 20:
        spacing = 5
    elif hrange < 120:
        spacing = 20
    elif hrange < 601:
        spacing = 60
    else:
        spacing = 300
    return spacing


def plot_twopoints(univariate, condition_label=None, trefs=[], ntrefs=4,
                   axe_xsize=6., axe_ysize=2.,
                   time_range=(None, None),
                   time_fractional_pad=.1,
                   counts_range=(None, None),
                   counts_fractional_pad=.1,
                   corr_range=(None, None),  # auto
                   corr_fractional_pad=.1,
                   delta_t_max=None,
                   show_exp_decay=None,
                   show_legend=True,
                   show_cdt_details_in_legend=False,
                   save=False, ext='.png'):
    """Plot two-point functions.

    Parameters
    ----------
    univariate : :class:`Univariate` instance
    condition_label : str (default None)
        must be the repr of a given FilterSet
    trefs : flist of floats
        indicate the times that you would like to have as references
        if left empty, reference times will be computed automatically
    ntrefs : int
        if trefs is empty, number of times of reference to display
    axe_xsize : float (default 6)
        size of the x-axis (inches)
    axe_ysize : float (default 2.)
        size if a single ax y-axis (inches)
    time_range : couple of floats (default (None, None))
        specifies (left, right) bounds
    time_fractional_pad : float (default .1)
        fraction of x-range to add as padding
    counts_range : couple of floats (default (None, None))
        specifies range for the Counts y-axis
    counts_fractional_pad : float (default .2)
        fractional amount of y-range to add as padding
    corr_range : couple of floats (default (None, None))
        sepcifies range for the Average y-axis
    corr_fractional_pad : couple of floats (default .2)
        fractional amounts of range to padding
    delta_t_max : float (default None)
        when given, bottom plot will be using this max range symmetrically;
        otherwise, will use the largest intervals found in data (often too
        large to see something)
    show_exp_decay : float (default None)
        when a floating point number is passed, a light exponential decay
        curve is plotted for each tref
     show_legend : bool {True, False}
        print out legend
    show_cdt_details_in_legend : bool {False, True}
        show details about filters
    save : bool {False, True}
        whether to save figure at canonical path
    ext : str {'.png', '.pdf'}
        extension to be used when saving figure
    """
    obs = univariate.obs
    # get priod from eval times
    if len(univariate.eval_times) > 0:
        period = univariate.eval_times[1] - univariate.eval_times[0]
    # or from experiment metadata
    else:
        period = univariate.exp.period
    fig, axs = plt.subplots(3, 1, figsize=(axe_xsize, 3*axe_ysize))

    # define time label
    if univariate.obs.mode != 'dynamics' and univariate.obs.timing == 'g':
        timelabel = 'Generations'
        if univariate.obs.tref is not None:
            timelabel += ' (since tref {})'.format(univariate.obs.tref)
    else:
        timelabel = 'Time (minutes)'

    # choice of index/indices for time of reference
    times = univariate['master'].time
    npoints = len(times)
    if not trefs:
        logging.info('Determining trefs...')
        di = npoints // ntrefs + 1
        indices = np.arange(0, npoints, di, dtype=int)
        trefs = times[indices]
        logging.info(trefs)

    ax01_mins = []
    ax01_maxs = []
    
    all_counts = []
    all_corr = []

    default_lw = mpl.rcParams['lines.linewidth']
    handles = []

    conditions = ['master', ] + univariate.cset
    for index, cdt in enumerate(conditions):
        if cdt == 'master':
            c_repr = 'master'
            c_label = 'all samples'
            lw = default_lw + 1
            lt = '-'
            alpha = .8
        elif cdt.label == condition_label or str(cdt) == condition_label or repr(cdt) == condition_label:
            c_repr = repr(cdt)
            if show_cdt_details_in_legend:
                c_label = str(cdt)
            else:
                c_label = cdt.label
            lw = default_lw
            lt = '--'
            alpha = .6
        # we plot master and one condition if given, not more...
        else:
            continue

        times = univariate[c_repr].time
        counts = univariate[c_repr].count_two
        corr = univariate[c_repr].autocorr
        var = np.diagonal(corr)

        valid = counts != 0
        
        latex_ref = '{{\mathrm{{ref}}}}'
        if obs.timing == 'g':
            prefix = 'g'
            units = ''
        else:
            prefix = 't'
            units = 'mins'

        for tref in trefs:
            # this tref may not be in conditioned data (who knows)
            if np.amin(np.abs(times - tref)) > period:
                continue
            index = np.argmin(np.abs(times - tref))
            if obs.timing == 'g':
                lab = '{:d}'.format(tref)
            else:
                lab = '{:.0f}'.format(tref)
            line_label = r'$ {}_{} = {}$ {} ({})'.format(prefix, latex_ref, lab, units, c_label)

            ok = np.where(counts[index, :] > 0)
#            if len(ok[0]) == 0:
#                continue
            # time limits
            xmin, xmax = np.nanmin(times[ok]), np.nanmax(times[ok])
            ax01_mins.append(xmin)
            ax01_maxs.append(xmax)
            dat, = axs[0].plot(times[ok], counts[index, :][ok],
                      ls=lt, lw=lw, alpha=alpha, label=line_label)
            handles.append(dat)
            all_counts.extend(counts[index, :][ok])
            color = dat.get_color()
            axs[0].plot((tref, tref), (0, counts[index, index]),
                        ls=':', color=color)

            axs[1].axhline(0, ls='-', color='C7', alpha=.3)  # thin line at 0
            dat, = axs[1].plot(times[valid[index, :]],
                           corr[index, :][valid[index, :]]/var[index],
                           ls=lt, lw=lw, alpha=alpha)
            all_corr.extend(corr[index, :][valid[index, :]]/var[index])
            color = dat.get_color()
            
            axs[1].axvline(tref, ymin=0.1, ymax=0.9, ls=':', color=color)
            
            axs[2].axhline(0, ls='-', color='C7', alpha=.3)  # thin line at 0
            axs[2].plot(times[valid[index, :]] - tref,
                    corr[index, :][valid[index, :]]/var[index], ls=lt, lw=lw, alpha=alpha)
    
    # xaxis limits
    if ax01_mins:
        left = np.amin(ax01_mins)
    else:
        left = None
    if ax01_maxs:
        right = np.amax(ax01_maxs)
    else:
        right = None
    if time_range[0] is not None:
        left = time_range[0]
    if time_range[1] is not None:
        right = time_range[1]
    hrange = right - left
    for ax in axs[:2]:
        ax.set_xlim(left=left - time_fractional_pad*hrange,
                    right=right + time_fractional_pad*hrange)
    # bottom plot : try to zoom over provided range
    if delta_t_max is not None:
        axs[2].set_xlim(left=-delta_t_max, right=delta_t_max)
    # if not provided, compute automatic ranges (not pretty usually)
    elif left is not None and right is not None:
        axs[2].set_xlim(left=-hrange, right=hrange)

    # add exponential decay
    if show_exp_decay is not None:
        tt = np.linspace(left, right, 100)
        dd = np.linspace(-hrange, hrange, 100)
        lab = r'$t_{{\mathrm{{decay}}}} = {:.1f}$ {}'.format(1./show_exp_decay, units)
        for tref in trefs:
            axs[1].plot(tt, np.exp(-show_exp_decay * np.abs(tt - tref)),
                        ls='-.', color='C7', alpha=.7)
        dec, = axs[2].plot(dd, np.exp(-show_exp_decay * np.abs(dd)),
                    ls='-.', color='C7', alpha=.5, label=lab)
        all_corr.extend(np.exp(-show_exp_decay * np.abs(dd)))
        handles.append(dec)

    # yaxis limits
    max_count = np.nanmax(all_counts)
    axs[0].set_ylim(bottom=0, top=max_count*(1. + counts_fractional_pad))  # counts

    # average
    bottom, top = np.nanmin(all_corr), np.nanmax(all_corr)
    if corr_range[0] is not None:
        bottom = corr_range[0]
    if corr_range[1] is not None:
        top = corr_range[1]
    vrange = top - bottom
    for ax in axs[1:]:
        ax.set_ylim(bottom=bottom - corr_fractional_pad*vrange,
                        top=top + corr_fractional_pad*vrange)

    labels = [h.get_label() for h in handles]
    axs[-1].legend(handles=handles, labels=labels, loc='upper left',
                   bbox_to_anchor=(0, -.6/axe_ysize), labelspacing=0.2)  # reduce labelspacing because of LaTeX

    # ticks and labels
    for ax in axs:
        ax.xaxis.set_major_locator(ticker.MultipleLocator(_spacing_trange(hrange)))
    # counts
    axs[0].yaxis.set_major_locator(MaxNLocator(nbins=3, integer=True))
    for ax in axs[1:]:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=2))
    for ax in axs[1:]:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=3))
    formatter = ScalarFormatter(useMathText=True, useOffset=False)
    formatter.set_powerlimits((-2, 4))
    for ax in axs:
        ax.yaxis.set_major_formatter(formatter)
        t = ax.yaxis.get_offset_text()
        plt.draw()
        msg = t.get_text()
        ax.text(0, .95, msg, ha='left', va='top', transform=ax.transAxes)
        t.set_visible(False)

    axs[0].tick_params(axis='x', direction='in', bottom='on', labelbottom='on', pad=-10)
    axs[1].tick_params(axis='x', direction='in', bottom='on', labelbottom='on', pad=-10)
    axs[2].set_xlabel(timelabel, x=.95, horizontalalignment='right',
                      fontsize='large')
    
    # hide intermediate x axis
    for ax in axs[:1]:
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='x', colors='C7')
    for ax in axs[1:2]:
        ax.spines['top'].set_color('C7')

    # ylabels
    axs[0].set_ylabel(r'# $\langle t_{\mathrm{ref}} | t \rangle$',
                      fontsize='large')
    axs[1].set_ylabel(r'$a(t_{\mathrm{ref}}, t)$',
                      fontsize='large')
    axs[2].set_ylabel(r'$a(t_{\mathrm{ref}}, t- t_{\mathrm{ref}})$',
                      fontsize='large')

    # title
    latex_obs = obs.as_latex_string
    axs[0].text(0.5, 1+.2/axe_ysize,
                r' Autocorrelation: {}'.format(latex_obs),
                size='x-large',
                horizontalalignment='center',
                verticalalignment='bottom',
                transform=axs[0].transAxes)
    fig.subplots_adjust(hspace=0.)
    # save fig at canonical path
    if save:
        # export data files if not existing yet
        try:
            obs_path = univariate._get_obs_path(write=False)
        except text.MissingFolderError:
            univariate.write_text()
        if condition_label is None:
            univc = univariate.master
        else:
            univc = univariate[condition_label]
        cdt_path = univc._get_path()
        bname = 'plot_twopoints_' + univariate.obs.name + '_' + univariate.region.name + ext
        fname = os.path.join(cdt_path, bname)
        fig.savefig(fname, bbox_inches='tight', pad_inches=0)
    return fig


def plot_stationary(stationary, show_cdts='all',
                    axe_xsize=6., axe_ysize=2.,
                    time_range=(None, None),
                    time_fractional_pad=.1,
                    counts_range=(None, None),
                    counts_fractional_pad=.1,
                    corr_range=(None, None),  # auto
                    corr_fractional_pad=.1,
                    fitlog=False, epsilon=0.1, fitting_time_max=None,
                    ref_decay=None,
                    interval_max=None, save=False, ext='.png'):
    """Plot stationary autocorrelation.

    Parameters
    ----------
    stationary : StationaryUnivariate or StationaryBivariate instance
    fitlog : bool {False, True}
        whether to fit initial decay with an exponential decay
    epsilon : float (default 0.1)
        threshold to ensure 'close to zero' in the fitting procedure, tries
        to perform the fit on autocorrelation values larger than epsilon.
        Default value makes the fit over one decade.
    fitting_time_max : float (default None)
        maximum time range to make fit when provided. Fit will be performed
        over the smaller region defined by both epsilon and fitting_time_max
    ref_decay : float (default None)
        whether to plot an exponential decay with corresponding rate
        exp(-rate * t)
    interval_max : float (default None)
        fixes the largest interval (instead of largest interval found in data)
    save : bool {False, True}
        whether to save plot at canonical path
    ext : str {'.png', '.pdf'}
        extension used for file

    Returns
    -------
    fig : Figure instance
    """
    if not (isinstance(stationary, StationaryUnivariate) or
            isinstance(stationary, StationaryBivariate)):
        msg = ('Input is not an instance of '
               '{}'.format(StationaryUnivariate) + 'or of '
               '{}'.format(StationaryBivariate))
        raise TypeError(msg)
    if isinstance(stationary, StationaryUnivariate):
        obs = stationary.obs
    elif isinstance(stationary, StationaryBivariate):
        obs = [uni.obs for uni in stationary.univariates]
#    if fitlog and isinstance(stationary, StationaryUnivariate):
#        nplots = 3
#    else:
#        nplots = 2
    nplots = 3
    fig = plt.figure(figsize=(axe_xsize, (nplots + 1)*axe_ysize))
    gs = gridspec.GridSpec(nplots + 1, 1)
    ax1 = fig.add_subplot(gs[0])
    if nplots == 2:
        ax2 = fig.add_subplot(gs[1:])
    elif nplots == 3:
        ax2 = fig.add_subplot(gs[1:-1])
        ax3 = fig.add_subplot(gs[-1:])

    additional_handles = []

    # build condition list
    conditions = ['master', ]  # list of conditions to be plotted
    if show_cdts == 'all':
        conditions = ['master', ] + stationary.cset
    elif show_cdts == 'master':
        pass
    elif isinstance(show_cdts, collections.Iterable):
        for item in show_cdts:
            _append_cdt(stationary.univariate, item, conditions)
    else:
        _append_cdt(stationary.univariate, show_cdts, conditions)

    tleft = np.infty
    tright = - np.infty
    xleft = np.infty
    xright = - np.infty
    for index, condition in enumerate(conditions):
        letter = string.ascii_lowercase[index]

        if condition == 'master':
            cdt_repr = 'master'
            cdt_str = 'all samples'
        else:
            cdt_repr = repr(condition)
            cdt_str = str(condition)
        array = stationary[cdt_repr].array
        nonzero = np.where(array['counts'] > 1)  # 1 sample does not have std
        dts = array['time_interval'][nonzero]
        if np.nanmax(dts) > tright:
            tright = np.nanmax(dts)
        if np.nanmin(dts) < tleft:
            tleft = np.nanmin(dts)
        counts = array['counts'][nonzero]

        if isinstance(stationary, StationaryUnivariate):
            corr = array['auto_correlation'][nonzero]
        else:
            corr = array['cross_correlation'][nonzero]
        try:
            dev = array['std_dev'][nonzero]
        except ValueError:
            dev = None

        # counts
        label = '({}) {}'.format(letter, cdt_str)
        line, = ax1.plot(dts, counts, label=label)
        col = line.get_color()  # usefule for later stage

        # autocorrelation: divide by variance
        if isinstance(stationary, StationaryUnivariate):
            norm = corr[0]
        # cross-correlation: divide covariance by product of standard devs
        elif isinstance(stationary, StationaryBivariate):
            prod = 1.
            for single in stationary.univariates:
                prod *= np.sqrt(single[cdt_repr].stationary.autocorr[0])
            norm = prod
        dat, = ax2.plot(dts, corr/norm, color=col,
                        label=label)
        if dev is not None:
            se = 2.58 * dev / np.sqrt(counts)
            ci = ax2.fill_between(dts, (corr-se)/norm, (corr+se)/norm,
                                  facecolor=col, alpha=.5,
                                  label='.99 C.I. for ({})'.format(letter))
            additional_handles.append(ci)

        ax2.axhline(0, ls='--', color='k', alpha=.5)
        ax2.axvline(0, ls=':', color='k', alpha=.5)

        # plot exponential fit for autocorrelation
        if fitlog and isinstance(stationary, StationaryUnivariate):
            # go though positive values larger than epsilon
            for arg, val in enumerate(corr[:]/corr[0]):
                if val < epsilon:
                    break
            argmax = arg
            if fitting_time_max is not None:
                argmax = np.nanargmin(np.abs(dts - fitting_time_max))
            # choose minimum between arg and argmax
            if argmax > arg:
                argmax = arg
            if argmax < 2:
                warnings.warn('Not enough points to perform fitlog')
            else:
                xdata = dts[:argmax]
                if np.nanmax(xdata) > xright:
                    xright = np.nanmax(xdata)
                if np.nanmin(xdata) < xleft:
                    xleft = np.nanmin(xdata)
                ydata = corr[:argmax]/corr[0]
                ci = se[:argmax]/norm
                slope, intercept = np.polyfit(xdata, np.log(ydata), 1)
                tt = np.linspace(0, xdata[argmax-1])
                ax3.plot(xdata, ydata, ls='none', marker='o', ms=8,
                        markerfacecolor='none',
                        markeredgecolor=col)
                ax3.fill_between(xdata, ydata-ci, ydata+ci, facecolor=col,
                                 alpha=.5)
                ax3.plot(tt, np.exp(slope * tt + intercept), color=col,
                         label=r'$\tau = {:.1f}$ mins'.format(-1./slope))

    if ref_decay is not None:
        tt = np.linspace(tleft, tright, 50)
        yy = np.exp(-ref_decay*np.abs(tt))
        ref, = ax2.plot(tt, yy, '--', color='k', alpha=.5,
                        label=r'$\tau = {:.1f}$ mins'.format(1./ref_decay))
        additional_handles.append(ref)
        if fitlog:
            tt = np.linspace(xleft, xright, 20)
            yy = np.exp(-ref_decay*tt)
            ax3.plot(tt, yy, '--', color='k', alpha=.5,
                     label=r'$\tau = {:.1f}$ mins'.format(1./ref_decay))
    
    if interval_max is not None and interval_max < tright:
        _right = interval_max
        if isinstance(stationary, StationaryUnivariate):
            _left = 0.
        else:
            _left = -interval_max
        ax1.set_xlim(left=_left, right=_right)
        ax2.set_xlim(left=_left, right=_right)
        if fitlog:
            ax3.set_xlim(left=_left, right=_right)

    formatter = ScalarFormatter(useMathText=True, useOffset=False)
    formatter.set_powerlimits((-2, 3))
    ax1.yaxis.set_major_formatter(formatter)
    # ylabels
    ax1.set_ylabel(r'Samples', fontsize='large')
    ax2.set_ylabel(r'Autocorr. $g_s(\Delta t)$', fontsize='large')

    ax2.set_ylim(bottom=-1.2, top=1.2)
    if len(conditions) > 2:
        ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    else:
        ax1.legend(loc=0)
    if nplots == 2:
        ax = ax2
    elif nplots == 3:
        ax = ax3
    if isinstance(stationary, StationaryUnivariate):
        ob = obs
    elif isinstance(stationary, StationaryBivariate):
        ob = obs[0]
    if ob.timing == 'g':
        units = 'generations'
    else:
        units = 'mins'
    ax.set_xlabel(r'$\Delta t$ ({})'.format(units), fontsize='large')

    # confidence interval legend
    if additional_handles:
        labels = [item.get_label() for item in additional_handles]
        if len(conditions) > 2:
            ax2.legend(handles=additional_handles, labels=labels,
                       bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        else:
            ax2.legend(handles=additional_handles, labels=labels, loc=0)

    if fitlog and isinstance(stationary, StationaryUnivariate):
        ax3.set_yscale('log')
        ax3.set_xlim(ax2.get_xlim())
        ax3.legend(loc=0)
    # writting observable
    # case: obs is a single observable
    if isinstance(stationary, StationaryUnivariate):
        latex_obs = obs.as_latex_string
        msg = 'autocorrelation'
    # case: obs is a couple of observables
    else:
        latex_obs = ', '.join([ob.as_latex_string for ob in obs])
        msg = 'covariance'

    ax1.text(0.5, 1.1, r' Stationary {} for {}'.format(msg, latex_obs),
             size='large',
             horizontalalignment='center',
             verticalalignment='bottom',
             transform=ax1.transAxes)
    if save:
        # get univariate instance to get path where to save figure
        bname = 'plot_stationary_'
        try:
            obs_path = stationary._get_obs_path(write=False)
        except text.MissingFolderError:
            stationary.write_text()
            obs_path = stationary._get_obs_path(write=False)
        obsname = os.path.basename(obs_path)
        bname += obsname + '_'
        bname += stationary.region.name + ext
        fname = os.path.join(obs_path, bname)
        fig.savefig(fname, bbox_inches='tight', pad_inches=0)
    return fig
