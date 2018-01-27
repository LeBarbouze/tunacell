#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module defines plotting functions for the statistics of the dynamics.
"""
from __future__ import print_function

import os
import numpy as np
import collections
import logging

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import ticker
import matplotlib.gridspec as gridspec

from tunacell.filters.main import FilterSet
from tunacell.stats.single import Univariate, StationaryUnivariate
from tunacell.stats.two import StationaryBivariate
from tunacell.io import text

from .helpers import _set_axis_limits, _set_timelabel, _set_time_axis_ticks


# few variables that will be used through all functions
default_fontsize = mpl.rcParams['font.size']
default_lw = mpl.rcParams['lines.linewidth']


def _set_condition_list(univariate, show_cdts='master'):
    """Set the list of conditions to show

    Parameters
    ----------
    show_cdts : str or FilterSet or iterable on these (default 'master')
        the conditions to plot, use 'all' for all conditions in univariate
    univariate : Univariate instance
        conditions will be matched against conditions stored in univariate

    Returns
    -------
    list of FilterSet (conditions) to show
    """
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
    return conditions


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
            # TODO : compare also cdt.label
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
                  use_obs_name=None,
                  save=False, user_path=None, ext='.png',
                  verbose=False):
    """Plot one point statistics: counts, average, abd variance.

    One point functions are plotted for each condition set up in *show_cdts*
    argument: 'all' for all conditions, or the string representation (or label)
    of a particuler condition (or a list thereof).

    Parameters
    ----------
    univariate : Univariate instance
    show_cdts : str (default 'all')
        must be either 'all', or 'master', or the repr of a condition, or a
        list thereof
    show_ci : bool {False, True}
        whether to show 99% confidence interval
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
    use_obs_name : str (default None)
        when filled, the plot title will use this observable name instead
        of looking for the observable registered name
    save : bool {False, True}
        whether to save plot
    user_path : str (default None)
        user defined path where to save figure; default is canonical path
        (encouraged)
    ext : str {'.png', '.pdf'}
        extension to be used when saving file
    verbose : bool {False, True}
    """
    if not isinstance(univariate, Univariate):
        raise TypeError('Input is not {}'.format(Univariate))

    fig, axs = plt.subplots(3, 1, figsize=(axe_xsize, 3*axe_ysize))

    obs = univariate.obs
    timelabel = _set_timelabel(obs)  # define time label

    main_handles = []  # main legend
    ci_handles = []  # additional legend (TODO: check if necessary)

    all_times = []
    all_counts = []
    all_average = []
    all_variance = []

    # build condition list
    conditions = _set_condition_list(univariate, show_cdts)

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
        se = 2.58 * std / np.sqrt(counts)  # standard error 99% CI Gaussian
#        var = np.diagonal(univariate[c_repr].autocorr)

        line_counts, = axs[0].plot(times, counts, alpha=alpha, lw=lw,
                               label='{}'.format(c_label))
        main_handles.append(line_counts)
        color = line_counts.get_color()

        average, = axs[1].plot(times, mean, color=color, alpha=0.8, lw=lw, label=c_label)
        if show_ci:
            fill_std = axs[1].fill_between(times, mean-se, mean+se,
                                           facecolor=color, alpha=alpha_fill)
            ci_handles.append(fill_std)
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

    # print vertical line at tref
    if obs.timing != 'g' and isinstance(obs.tref, float):
        for ax in axs:
            vtref = ax.axvline(univariate.obs.tref, color='C7', ls='--',
                                alpha=.5, label='reference time in obs')
        main_handles.append(vtref)  # only the last one

    # ## limits and ticks ##
    # xaxis
    for ax in axs:
        left, right = _set_axis_limits(ax, all_times, which='x', pad=time_fractional_pad,
                                       force_range=time_range)
    # locator
    locator = _set_time_axis_ticks(axs[0], obs, bounds=(left, right))
    for ax in axs:
        ax.xaxis.set_major_locator(locator)

    # yaxis limits
    _set_axis_limits(axs[0], all_counts, which='y', pad=counts_fractional_pad,
                     force_range=counts_range)
    axs[0].yaxis.set_major_locator(ticker.MaxNLocator(nbins=3, integer=True))
    # average
    _set_axis_limits(axs[1], all_average, which='y', pad=average_fractional_pad,
                     force_range=average_range)
    axs[1].yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    # variance
    _set_axis_limits(axs[2], all_variance, which='y', pad=variance_fractional_pad,
                     force_range=variance_range)
    axs[2].yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))
    # tick formatter
    formatter = ticker.ScalarFormatter(useMathText=True, useOffset=False)
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
                      fontsize='medium')

    # hide intermediate x axis
    for ax in axs[:2]:
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='x', colors='C7')
    for ax in axs[1:]:
        ax.spines['top'].set_color('C7')

    axs[0].set_ylabel('Counts', fontsize='medium')
    axs[1].set_ylabel('Average', fontsize='medium')
    axs[2].set_ylabel('Variance', fontsize='medium')

    # ## legend ##
    # C.I.
    if ci_handles:
        ci = ci_handles[0]
#        ci.set_color('C7')
        ci.set_label('.99 C.I.')
        main_handles.append(ci)

    handles = main_handles[:]
    labels = [h.get_label() for h in handles]
    if show_legend:
        axs[-1].legend(handles=handles, labels=labels, loc='upper left',
                       bbox_to_anchor=(0, -.5/axe_ysize))

    # title
    latex_obs = obs.latexify(use_name=use_obs_name)
    axs[0].text(0.5, 1+.2/axe_ysize,
                r'{}'.format(latex_obs),
                size='large',
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
        fig.savefig(fname, bbox_inches='tight', pad_inches=0)
        if verbose:
            print('Figure saved as {}'.format(fname))
    return fig


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
                   use_obs_name=None,
                   save=False, ext='.png', verbose=False):
    """Plot two-point functions: counts and autocorrelation functions.

    These plots are able to show only one extra condition with 'master', and
    are plotted for a set of time of references.

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
    use_obs_name : str (default None)
        when filled, the plot title will use this observable name instead
        of looking for the observable registered name
    save : bool {False, True}
        whether to save figure at canonical path
    ext : str {'.png', '.pdf'}
        extension to be used when saving figure
    verbose : bool {False, True}
    """
    obs = univariate.obs
    timelabel = _set_timelabel(obs)  # define time label
    # get priod from eval times
    if len(univariate.eval_times) > 0:
        period = univariate.eval_times[1] - univariate.eval_times[0]
    # or from experiment metadata
    else:
        period = univariate.exp.period
    fig, axs = plt.subplots(3, 1, figsize=(axe_xsize, 3*axe_ysize))

    # choice of index/indices for time of reference
    times = univariate['master'].time
    npoints = len(times)
    if not trefs:
        logging.info('Determining trefs...')
        di = npoints // ntrefs + 1
        indices = np.arange(0, npoints, di, dtype=int)
        trefs = times[indices]
        logging.info(trefs)

    all_times = []
    all_counts = []
    all_corr = []

    handles = []

    # prep work for latex printing
    latex_ref = '{{\mathrm{{ref}}}}'
    if obs.timing == 'g':
        prefix = 'g'
        units = ''
    else:
        prefix = 't'
        units = 'mins'

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
            all_times.extend(times[ok])
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

    # ## limits and ticks ##
    # xaxis
    for ax in axs[:2]:
        left, right = _set_axis_limits(ax, all_times, which='x',
                                       pad=time_fractional_pad,
                                       force_range=time_range)
        hrange = right - left
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    # bottom plot : try to zoom over provided range
    if delta_t_max is not None:
        axs[2].set_xlim(left=-delta_t_max, right=delta_t_max)
    # if not provided, compute automatic ranges (not pretty usually)
    else:
        axs[2].set_xlim(left=-hrange, right=hrange)
    axs[2].xaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    # add exponential decay
    if show_exp_decay is not None:
        tt = np.linspace(left, right, 100)
        dd = np.linspace(-hrange, hrange, 100)
        lab = r'$t_{{\mathrm{{decay}}}} = {:.1f}$ {}'.format(1./show_exp_decay, units)
        for tref in trefs:
            axs[1].plot(tt, np.exp(-show_exp_decay * np.abs(tt - tref)),
                        ls='-.', color='C7', alpha=.7)
        dec, = axs[2].plot(dd, np.exp(-show_exp_decay * np.abs(dd)),
                    ls='-.', color='C7', alpha=.7, label=lab)
        all_corr.extend(np.exp(-show_exp_decay * np.abs(dd)))
        handles.append(dec)

    # ## yaxis limits ##
    # counts
    _set_axis_limits(axs[0], all_counts, which='y', pad=counts_fractional_pad,
                     force_range=counts_range)
    axs[0].yaxis.set_major_locator(ticker.MaxNLocator(nbins=3, integer=True))

    # corr
    for ax in axs[1:]:
        _set_axis_limits(ax, all_corr, which='y', pad=corr_fractional_pad,
                         force_range=corr_range)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3))

    # legend
    labels = [h.get_label() for h in handles]
    axs[-1].legend(handles=handles, labels=labels, loc='upper left',
                   bbox_to_anchor=(0, -.5/axe_ysize), labelspacing=0.2)  # reduce labelspacing because of LaTeX

    formatter = ticker.ScalarFormatter(useMathText=True, useOffset=False)
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
                      fontsize='medium')

    # hide intermediate x axis
    for ax in axs[:1]:
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(axis='x', colors='C7')
    for ax in axs[1:2]:
        ax.spines['top'].set_color('C7')

    # ylabels
    axs[0].set_ylabel(r'# $\langle t_{\mathrm{ref}} | t \rangle$',
                      fontsize='medium')
    axs[1].set_ylabel(r'$a(t_{\mathrm{ref}}, t)$',
                      fontsize='medium')
    axs[2].set_ylabel(r'$a(t_{\mathrm{ref}}, t- t_{\mathrm{ref}})$',
                      fontsize='medium')

    # title
    latex_obs = obs.latexify(use_name=use_obs_name)
    axs[0].text(0.5, 1+.2/axe_ysize,
                r'{}'.format(latex_obs),
                size='large',
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
        bname = 'plot_twopoints_' + obs.name + '_' + univariate.region.name + ext
        fname = os.path.join(cdt_path, bname)
        fig.savefig(fname, bbox_inches='tight', pad_inches=0)
        if verbose:
            print('Figure saved as {}'.format(fname))
    return fig


def plot_stationary(stationary, show_cdts='all',
                    axe_xsize=6., axe_ysize=2.,
                    time_range=(None, None),
                    time_fractional_pad=.1,
                    time_guides=[0., ],
                    counts_range=(None, None),
                    counts_fractional_pad=.1,
                    corr_range=(None, None),  # auto
                    counts_logscale=False,
                    corr_fractional_pad=.1,
                    corr_logscale=False,
                    corr_guides=[0., ],
                    show_exp_decay=None,
                    show_legend=True, show_cdt_details_in_legend=False,
                    use_obs_name=None,
                    save=False, ext='.png', verbose=False):
    """Plot stationary autocorrelation.

    Parameters
    ----------
    stationary : StationaryUnivariate or StationaryBivariate instance
    axe_xsize : float (default 6)
        size (in inches) of the x-axis
    axe_ysize : float (default 2)
        size (in inches) of the individual y-axis
    time_range : couple of floats
        bounds for time (x-axis)
    time_fractional_pad : float
        fractional padding for x-axis
    counts_range : couple of ints
        bounds for counts axis
    counts_fractional_pad : float
        fractional padding for counts axis
    corr_range : couple of floats
        bounds for correlation values
    counts_logscale : bool {False, True}
        use logscale for counts axis
    corr_fractional_pad : float
        fractional padding for correlation values
    corr_logscale : bool {False, True}
        use logscale for correlation values (symlog is used to display
        symmetrically negative values)
    corr_guides : list of float
        values where to plot shaded grey horizontal lines
    show_exp_decay : float (default None)
        whether to plot an exponential decay with corresponding rate
        exp(-rate * t)
    save : bool {False, True}
        whether to save plot at canonical path
    use_obs_name : str (default None)
        when filled, the plot title will use this observable name instead
        of looking for the observable registered name
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
        timelabel = _set_timelabel(obs, use_tref=False)
    elif isinstance(stationary, StationaryBivariate):
        obs = [uni.obs for uni in stationary.univariates]
        timelabel = _set_timelabel(obs[0], use_tref=False)
    if 'minutes' in timelabel:
        units = 'mins'
        prefix = 't'
    else:
        units = ''  # generations are used
        prefix = 'g'
    timelabel = r'$\Delta$'+timelabel

    nplots = 2
    fig = plt.figure(figsize=(axe_xsize, (nplots + 1)*axe_ysize))
    gs = gridspec.GridSpec(nplots + 1, 1)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1:])

    # build condition list
    if isinstance(stationary, StationaryUnivariate):
        conditions = _set_condition_list(stationary.univariate, show_cdts=show_cdts)
    elif isinstance(stationary, StationaryBivariate):
        conditions = []
        conditions_0 = _set_condition_list(stationary.univariates[0], show_cdts=show_cdts)
        conditions_1 = _set_condition_list(stationary.univariates[1], show_cdts=show_cdts)
        # intersect
        for cdt in conditions_0:
            if cdt in conditions_1:
                conditions.append(cdt)

    all_times = []
    all_counts = []
    all_corrs = []

    main_handles = []  # for legend
    ci_handles = []

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

        array = stationary[c_repr].array
        nonzero = np.where(array['counts'] > 1)  # 1 sample does not have std
        dts = array['time_interval'][nonzero]
        all_times.extend(dts)
        counts = array['counts'][nonzero]
        all_counts.extend(counts)

        if isinstance(stationary, StationaryUnivariate):
            corr = array['auto_correlation'][nonzero]
        else:
            corr = array['cross_correlation'][nonzero]
        try:
            dev = array['std_dev'][nonzero]
        except ValueError:
            dev = None

        # counts
        label = '{}'.format(c_label)
        line, = ax1.plot(dts, counts, lw=lw, alpha=alpha, label=label)
        main_handles.append(line)
        col = line.get_color()  # usefule for later stage

        # autocorrelation: divide by variance
        if isinstance(stationary, StationaryUnivariate):
            norm = corr[0]
        # cross-correlation: divide covariance by product of standard devs
        elif isinstance(stationary, StationaryBivariate):
            prod = 1.
            for single in stationary.univariates:
                prod *= np.sqrt(single[c_repr].stationary.autocorr[0])
            norm = prod
        dat, = ax2.plot(dts, corr/norm, color=col,
                       lw=lw, alpha=alpha, label=label)
        all_corrs.extend(corr/norm)
        if dev is not None:
            se = 2.58 * dev / np.sqrt(counts)
            ci = ax2.fill_between(dts, (corr-se)/norm, (corr+se)/norm,
                                      facecolor=col, alpha=alpha_fill,
                                      label='.99 C.I.')
            ci_handles.append(ci)
            all_corrs.extend((corr-se)/norm)
            all_corrs.extend((corr+se)/norm)

    # vertical lines for timing
    for val in time_guides:
        ax2.axvline(val, ls=':', color='C7', alpha=.5)
    # horizontal lines for correlation ref
    for val in corr_guides:
        ax2.axhline(val, ls=':', color='C7', alpha=.5)
    # ## limits and ticks ##
    # xaxis
    for ax in [ax1, ax2]:
        left, right = _set_axis_limits(ax, all_times, which='x',
                                       pad=time_fractional_pad,
                                       force_range=time_range)
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    if show_exp_decay is not None:
        tt = np.linspace(left, right, 100)
        yy = np.exp(-show_exp_decay*np.abs(tt))
        lab = r'${}_{{\mathrm{{decay}}}} = {:.1f}$ {}'.format(prefix, 1./show_exp_decay, units)
        ref, = ax2.plot(tt, yy, '-.', color='C7', alpha=1,
                        label=lab)
        main_handles.append(ref)

    # ## yaxis limits ##
    # counts
    formatter = ticker.ScalarFormatter(useMathText=True, useOffset=False)
    formatter.set_powerlimits((-2, 4))
    if not counts_logscale:
        _set_axis_limits(ax1, all_counts, which='y', pad=counts_fractional_pad,
                         force_range=counts_range)
        ax1.yaxis.set_major_locator(ticker.MaxNLocator(nbins=3, integer=True))
        ax1.yaxis.set_major_formatter(formatter)
        t = ax.yaxis.get_offset_text()
        plt.draw()
        msg = t.get_text()
        ax1.text(0, .95, msg, ha='left', va='top', transform=ax.transAxes)
        t.set_visible(False)
    else:
        ax1.set_yscale('symlog', linthresh=1)

    # corr

    if not corr_logscale:
        bottom, top = _set_axis_limits(ax2, all_corrs, which='y',
                                       pad=corr_fractional_pad,
                                       force_range=corr_range)
        if top > 2 or bottom < -2:
            locator = ticker.MaxNLocator(nbins=5, integer=True)
        else:
            locator = ticker.FixedLocator([-1, -.5, 0., .5, 1])
        ax2.yaxis.set_major_locator(locator)
        ax2.yaxis.set_major_formatter(formatter)
        t = ax.yaxis.get_offset_text()
        plt.draw()
        msg = t.get_text()
        ax2.text(0, .95, msg, ha='left', va='top', transform=ax.transAxes)
        t.set_visible(False)
    else:
        ax2.set_yscale('symlog', linthreshy=0.1, linscaley=0.2,
                       subsy=[2, 3, 4, 5, 6, 7, 8, 9])
        if corr_range[0] is not None and corr_range[0] > 0.:
            ax2.set_ylim(bottom=corr_range[0])

    ax1.tick_params(axis='x', direction='in', bottom='on', labelbottom='on', pad=-10)
    ax2.set_xlabel(timelabel, x=.95, horizontalalignment='right',
                   fontsize='medium')

    # hide intermediate x axis
    ax1.spines['bottom'].set_visible(False)
    ax1.tick_params(axis='x', colors='C7')
    ax2.spines['top'].set_color('C7')

    # ylabels
    ax1.set_ylabel(r'Counts', fontsize='medium')
    if isinstance(stationary, StationaryUnivariate):
        ax2.set_ylabel(r'$\tilde{{a}}(\Delta {})$'.format(prefix), fontsize='medium')
    elif isinstance(stationary, StationaryBivariate):
        ax2.set_ylabel(r'$\tilde{{c}}(\Delta {})$'.format(prefix), fontsize='medium')

    # writting observable
    # case: obs is a single observable
    if isinstance(stationary, StationaryUnivariate):
        msg = '{}:{}'.format(obs.latexify(shorten_time_variable=True, use_name=use_obs_name),
                             obs.latexify(plus_delta=True, shorten_time_variable=True, use_name=use_obs_name))
    # case: obs is a couple of observables
    else:
        if use_obs_name is not None:
            if isinstance(use_obs_name, str):
                use_name_0 = use_obs_name
                use_name_1 = None
            else:
                if len(use_obs_name) == 1:
                    use_name_0 = use_obs_name[0]
                    use_name_1 = None
                else:
                    use_name_0 = use_obs_name[0]
                    use_name_1 = use_obs_name[1]
        else:
            use_name_0 = None
            use_name_1 = None
        msg = '{}:{}'.format(obs[0].latexify(shorten_time_variable=True,
                                             use_name=use_name_0),
                             obs[1].latexify(plus_delta=True, shorten_time_variable=True,
                                             use_name=use_name_1))

    ax1.text(0.5, 1+.2/axe_ysize, r'{}'.format(msg),
             size='large',
             horizontalalignment='center',
             verticalalignment='bottom',
             transform=ax1.transAxes)

    # ## legend ##
    # C.I.
    if ci_handles:
        ci = ci_handles[0]
        # ci.set_color('C7')
        ci.set_label('.99 C.I.')
        main_handles.append(ci)

    handles = main_handles[:]
    labels = [h.get_label() for h in handles]
    if show_legend:
        ax2.legend(handles=handles, labels=labels, loc='upper left',
                       bbox_to_anchor=(0, -.25/axe_ysize), labelspacing=.2)

    fig.subplots_adjust(hspace=0)
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
        if verbose:
            print('Figure saved as {}'.format(fname))
    return fig
