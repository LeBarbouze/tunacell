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

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, ScalarFormatter
import matplotlib.transforms as transforms
import matplotlib.lines as mlines
import matplotlib.gridspec as gridspec

from tuna.filters.main import FilterSet
from tuna.stats.single import Univariate, StationaryUnivariate
from tuna.stats.two import StationaryBivariate
from tuna.io import text


class UnivariatePlot(object):
    """Class to make onepoint and twopoint plotting

    .. note:: Deprecated in tuna 0.0.7

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


def plot_onepoint(univariate, show_cdts='all', left=None, right=None,
                  mean_ref=None, mean_lims=(None, None),
                  show_ci=False,
                  var_ref=None, var_lims=(None, None),
                  save=False, user_path=None, ext='.png'):
    """Plot one point statistics.

    Parameters
    ----------
    univariate : Univariate instance
    show_cdts : str (default 'all')
        must be either 'all', or 'master', or the repr of a condition, or a
        list thereof
    left : float (default None)
        leftest value to display
    right : float (default None)
        rightest value to display
    mean_ref : float
        reference mean value: what user expect to see as sample average to
        compare with data
    mean_lims : (float, float) (default (None, None))
        required limits for sample mean axes.yaxis
    var_ref : float
        reference variance value: what user expect to see as sample variance to
        compare with data
    var_lims : (float, float) (default (None, None))
        required limits for sample variance axes.yaxis
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
    vsize = 2.5
    fig, axs = plt.subplots(3, 1, figsize=(6, 3*vsize))

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

    additional_handles = []

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

    for index, cdt in enumerate(conditions):
        letter = string.lowercase[index]
        if cdt == 'master':
            c_repr = 'master'
            c_label = 'all samples'
        else:
            c_repr = repr(cdt)
            c_label = str(cdt)

        ok = np.where(univariate[c_repr].count_one > 0)

        times = univariate[c_repr].time[ok]
        counts = univariate[c_repr].count_one[ok]
        mean = univariate[c_repr].average[ok]
        var = univariate[c_repr].var[ok]
        std = univariate[c_repr].std[ok]
        se = 2.58 * std / np.sqrt(counts)
#        var = np.diagonal(univariate[c_repr].autocorr)

        ax = axs[0]
        line_counts, = ax.plot(times, counts, alpha=0.8,
                               label='({}) {}'.format(letter, c_label))

        color = line_counts.get_color()

        ax = axs[1]
        ax.plot(times, mean, color=color, alpha=0.8, label=c_label)
        if show_ci:
            fill_std = ax.fill_between(times, mean-se, mean+se,
                                       facecolor=color, alpha=.3,
                                       label='.99 C.I. for ({})'.format(letter))
            # add onbly if empty (no need to repeat on conditions)
            additional_handles.append(fill_std)

        ax = axs[2]
        ax.plot(times, var, color=color, alpha=0.8, label=c_label)

        ax.set_ylim(bottom=0)

    # adding reference lines
    if mean_ref is not None:
        ax = axs[1]
        ax.axhline(mean_ref, ls='--', color='C7', alpha=.7)
    if var_ref is not None:
        ax = axs[2]
        ax.axhline(var_ref, ls='--', color='C7', alpha=.7)

    # xaxis limits
    for ax in axs:
        if left is not None:
            ax.set_xlim(left=left)
        if right is not None:
            ax.set_xlim(right=right)
    
    # yaxis limits
    axs[0].set_ylim(bottom=0)  # counts

    mean_min, mean_max = mean_lims
    if mean_min is not None or mean_max is not None:
        ax = axs[1]
        if mean_min is not None:
            ax.set_ylim(bottom=mean_min)
        if mean_max is not None:
            ax.set_ylim(top=mean_max)
    var_min, var_max = var_lims
    if var_min is not None or var_max is not None:
        ax = axs[2]
        if var_min is not None:
            ax.set_ylim(bottom=var_min)
        if var_max is not None:
            ax.set_ylim(top=var_max)

    # print vertical line at tref
    if univariate.obs.timing != 'g' and isinstance(univariate.obs.tref, float):
        for ax in axs:
            ax.axvline(univariate.obs.tref, color='C7', ls='--', alpha=.5)

    # ticks and labels
    # first axes locators are integers
    axs[0].yaxis.set_major_locator(MaxNLocator(integer=True))
    formatter = ScalarFormatter(useMathText=True, useOffset=False)
    formatter.set_powerlimits((-2, 4))
    for ax in axs:
        ax.yaxis.set_major_formatter(formatter)
        t = ax.yaxis.get_offset_text()
        plt.draw()
        msg = t.get_text()
        ax.text(0, .95, msg, ha='left', va='top', transform=ax.transAxes)
        t.set_visible(False)
    axs[0].tick_params(axis='x', direction='out', top='on',
                       labeltop='on')
    axs[0].tick_params(axis='x', direction='in', bottom='on',
                       labelbottom='off')
    axs[1].tick_params(axis='x', direction='in')
    axs[1].set_xticklabels([])
    axs[2].set_xlabel(timelabel, x=.95, horizontalalignment='right',
                      fontsize='large')
    axs[0].xaxis.set_label_position('top')
    axs[0].set_xlabel(timelabel, x=.95, horizontalalignment='right',
                      fontsize='large')

    axs[0].set_ylabel('Counts', fontsize='large')
    axs[1].set_ylabel('Average', fontsize='large')
    axs[2].set_ylabel('Variance', fontsize='large')
    # legend
    if len(conditions) > 2:
        axs[0].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    else:
        axs[0].legend(loc=0)
    if additional_handles:
        labels = [item.get_label() for item in additional_handles]
        if len(conditions) > 2:
            axs[1].legend(handles=additional_handles, labels=labels,
                          bbox_to_anchor=(1.05, 1), loc=2)
        else:
            axs[1].legend(handles=additional_handles, labels=labels, loc=0)

    # add info above first ax
#    info = '(binning interval {})'.format(univariate.indexify.binsize)
#    axs[2].text(0., -.2, info, va='top', transform=axs[2].transAxes)

    # title
    latex_obs = obs.as_latex_string
    axs[0].text(0.5, 1.3, r' One-point stats for {}'.format(latex_obs),
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
        bname = 'plot_onepoint_' + univ.region.name + ext
        fname = os.path.join(obs_path, bname)
        fig.savefig(fname, bbox_to_inches='tight', pad_inches=0)
    return fig


def plot_twopoints(univariate, condition_label=None, trefs=[], ntrefs=4,
                   trange=(-100., 100.), show_exp_decay=None,
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
    trange : couple of floats
        limits of x-axis (time)
    show_exp_decay : float (default None)
        when a floating point number is passed, a light exponential decay
        curve is plotted for each tref
    save : bool {False, True}
        whether to save figure at canonical path
    ext : str {'.png', '.pdf'}
        extension to be used when saving figure
    """
    obs = univariate.obs
    fig, axs = plt.subplots(3, 1, figsize=(6, 9))

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
        print('Determining trefs...')
        di = npoints // ntrefs + 1
        indices = np.arange(0, npoints, di, dtype=int)
        trefs = times[indices]
        print(trefs)

    ax01_mins = []
    ax01_maxs = []

    if condition_label is None or condition_label == 'master':
        labels = ['master', ]
    else:
        labels = ['master', condition_label]

    for c_label in labels:
        if c_label is None:
            continue
        if c_label == 'master':
            lt = '-'
        else:
            lt = '--'

        times = univariate[c_label].time
        counts = univariate[c_label].count_two
        corr = univariate[c_label].autocorr
        var = np.diagonal(corr)

        valid = counts != 0
        
        tref_latex_label = 't_{{\mathrm{{ref}}}}'
        gref_latex_label = 'g_{{\mathrm{{ref}}}}'

        for tref in trefs:
            # this tref may not be in conditioned data (who knows)
            if np.amin(np.abs(times - tref)) > 1.:
                continue
            index = np.argmin(np.abs(times - tref))
            if obs.timing == 'g':
                lab = '{:d}'.format(tref)
                line_label = r'$' + gref_latex_label + '=' + lab +'$'
            else:
                lab = '{:.0f} mins'.format(tref)
                line_label = r'$' + tref_latex_label + '=' + lab + '$'

            ax = axs[0]
            ok = np.where(counts[index, :] > 0)
            if len(ok[0]) == 0:
                continue
            # time limits
            xmin, xmax = np.nanmin(times[ok]), np.nanmax(times[ok])
            ax01_mins.append(xmin)
            ax01_maxs.append(xmax)
            dat, = ax.plot(times[ok], counts[index, :][ok], ls=lt, label=line_label)
            color = dat.get_color()
            ax.plot((tref, tref), (0, counts[index, index]),
                    ls=':', color=color)

            ax = axs[1]
            dat, = ax.plot(times[valid[index, :]],
                           corr[index, :][valid[index, :]]/var[index],
                           ls=lt)
            color = dat.get_color()
            
            ax.axvline(tref, ymin=0.1, ymax=0.9, ls=':', color=color)
            ax.axhline(0, ls='--', color='k')
            
            # add text to point to tref
            # figure out where
            if obs.timing == 'g':
                pos = times[index] + 0.1
            elif index < len(times) - 1:
                pos = times[index + 1]
            else:
                pos = times[index - 2]
            for ax in axs[:2]:
                # define heterogeneous transform
                # the x coords of this transformation are data, and the
                # y coord are axes
                trans = transforms.blended_transform_factory(ax.transData,
                                                             ax.transAxes)
                ax.text(pos, 0.05, lab, color=color, transform=trans)
            # xmin, xmax = ax.xaxis.get_data_interval()

            ax = axs[2]
            ax.plot(times[valid[index, :]] - tref,
                    corr[index, :][valid[index, :]]/var[index], ls=lt)
            ax.axhline(0, ls='--', color='k')
    
    # axes limits
    if ax01_mins:
        left = np.amin(ax01_mins)
    else:
        left = None
    if ax01_maxs:
        right = np.amax(ax01_maxs)
    else:
        right = None
    for ax in axs[:2]:
        ax.set_xlim(left=left, right=right)
    delta_left, delta_right = left-right, right-left
    axs[2].set_xlim(left=delta_left, right=delta_right)

    axs[0].set_ylim(bottom=0)  # counts
    
    # add exponential decay
    if show_exp_decay is not None:
        tt = np.linspace(left, right, 100)
        dd = np.linspace(delta_left, delta_right, 100)
        for tref in trefs:
            axs[1].plot(tt, np.exp(-show_exp_decay * np.abs(tt - tref)),
                        ls='-.', color='C7', alpha=.7)
        axs[2].plot(dd, np.exp(-show_exp_decay * np.abs(dd)),
                    ls='-.', color='C7', alpha=.7)

    # legends
    axs[0].legend(loc=2)

    master_line = mlines.Line2D([], [], color='C7', ls='-',
                                label='master (color: see top)')
    handles = [master_line, ]
    if condition_label is not None and condition_label != 'master':
        c_line = mlines.Line2D([], [], color='C7', ls='--',
                               label=condition_label)
        handles.append(c_line)
    if show_exp_decay is not None:
        exp_line = mlines.Line2D([], [], color='C7', ls='-.', label='exp decay')
        handles.append(exp_line)
    axs[2].legend(handles=handles, loc=2)

    # ticks and labels
    # first axes locators are integers
    axs[0].yaxis.set_major_locator(MaxNLocator(integer=True))
    axs[0].tick_params(axis='x', direction='out', top='on',
                       labeltop='on')
    axs[0].tick_params(axis='x', direction='in', bottom='on',
                       labelbottom='off')
    axs[1].tick_params(axis='x', direction='in')

    axs[2].set_xlabel('Delta' + timelabel, x=.95, horizontalalignment='right',
                      fontsize='large')
    axs[0].xaxis.set_label_position('top')
    axs[0].set_xlabel(timelabel, x=.95, horizontalalignment='right',
                      fontsize='large')

    # ylabels
    axs[0].set_ylabel(r'Samples $\langle t_{\mathrm{ref}} | t \rangle$',
                      fontsize='large')
    axs[1].set_ylabel(r'Autocorr. $g(t_{\mathrm{ref}}, t)$',
                      fontsize='large')
    axs[2].set_ylabel(r'Shifted $g(t_{\mathrm{ref}}, t- t_{\mathrm{ref}})$',
                      fontsize='large')

    latex_obs = obs.as_latex_string
    axs[0].text(0.5, 1.3, r' Autocorrelation fcts for {}'.format(latex_obs),
                size='large',
                horizontalalignment='center',
                verticalalignment='bottom',
                transform=axs[0].transAxes)
    fig.subplots_adjust(hspace=.1)
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
        bname = 'plot_twopoints_' + univariate.region.name + ext
        fname = os.path.join(cdt_path, bname)
        fig.savefig(fname, bbox_inches='tight', pad_inches=0)
    return fig


def plot_stationary(stationary, show_cdts='all',
                    fitlog=False, epsilon=0.1, ref_decay=None,
                    interval_max=None, save=False, ext='.png'):
    """Plot stationary autocorrelation.

    Parameters
    ----------
    stationary : StationaryUnivariate or StationaryBivariate instance
    fitlog : bool {False, True}
        whether to fit initial decay with an exponential decay
    epsilon : float
        threshold to ensure 'close to zero' in the fitting procedure
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
    if fitlog and isinstance(stationary, StationaryUnivariate):
        nplots = 3
    else:
        nplots = 2

    fig = plt.figure(figsize=(6, 2 * (nplots + 1)))
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
        letter = string.lowercase[index]

        if condition == 'master':
            cdt_repr = 'master'
            cdt_str = 'all samples'
        else:
            cdt_repr = repr(condition)
            cdt_str = str(condition)
        array = stationary[cdt_repr].array
        nonzero = np.where(array['count'] > 1)  # 1 sample does not have std
        dts = array['time_interval'][nonzero]
        if np.nanmax(dts) > tright:
            tright = np.nanmax(dts)
        if np.nanmin(dts) < tleft:
            tleft = np.nanmin(dts)
        counts = array['count'][nonzero]

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
        yy = np.exp(-ref_decay*tt)
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
        if isinstance(stationary, StationaryUnivariate):
            bname = 'stationary_'
        elif isinstance(stationary, StationaryBivariate):
            bname = 'stationary_crosscorrelation_'
        try:
            obs_path = stationary._get_obs_path(write=False)
        except text.MissingFolderError:
            stationary.write_text()
            obs_path = stationary._get_obs_path(write=False)
        bname += stationary.region.name + ext
        fname = os.path.join(obs_path, bname)
        fig.savefig(fname, bbox_inches='tight', pad_inches=0)
    return fig
