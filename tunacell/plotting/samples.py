#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module defines plotting functions for samples (colonies, lineages)
"""
from __future__ import print_function

import os
import warnings
import re
import datetime
import logging

import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.ticker as ticker

import numpy as np

from tunacell.base.observable import set_observable_list

from tunacell.plotting.timeseries import add_timeseries, add_data_statistics
from .helpers import _set_time_axis_ticks, _unroll_samples

try:
    from StringIO import StringIO  # python2
except ImportError:
    from io import StringIO  # python3


class SamplePlot(object):
    """General class to store Colony plots.

    Parameters
    ----------
    samples : {(list of) :class:`Colony` instance(s),
               (list of) :class:`Lineage` instance(s)}
        timeseries from these samples will be plotted
    parser : :class:`Parser` instance
    conditions : list of :class:`FilterSet` instances
    label : str (default None)
        title for file saving; when None, current date and time are used
    """

    def __init__(self, samples, parser=None, conditions=[], label=None):
        self.parser = parser
        self.cset = conditions
        all_filters = [parser.experiment.fset, ] + conditions
        self._all_filters = all_filters  # used later to get all raw, func obs
        self._input_samples = list(samples)
        self._samples = []  # will be filled at make_plot call
        self.obs = None  # will be filled at make_plot call
        today = datetime.datetime.today()
        if label is not None:
            self.label = label
        else:
            self.label = 's{}'.format(today.strftime('%Y%m%d%H%M%S'))

        self.fig = None  # Figure instance

        return

    def make_plot(self, obs, **kwargs):
        """Produce Figure and save plotted data.

        Parameters
        ----------
        obs : :class:`Observable` isntance
            observable to plot
        kwargs : keyword arguments
            check :fun:`plot_samples` doctring for valid keyword arguments
        """
        self.obs = obs
        raw_obs, func_obs = set_observable_list(obs, filters=self._all_filters)
        raw_obs = raw_obs
        func_obs = func_obs
        samples = []
        for lin in _unroll_samples(self._input_samples):
            ts = lin.get_timeseries(obs, raw_obs, func_obs, self.cset)
            samples.append((lin, ts))
        self._samples = samples
        fig = plot_samples(samples, obs,
                           conditions=self.cset,
                           parser=self.parser,
                           **kwargs)
        self.fig = fig
        return

    def data_as_text(self, sep='\t', cell_sep='\n',
                     timeseries_sep='\n\n', print_labels=True):
        """Produces text output for plotted TimeSeries.

        Parameters
        ----------
        sep : str (default '\t')
            column separator
        cell_sep : str (default '\n')
            how to separate cell chunks in each timeseries
            (default is one blank line)
        timeseries_sep : str (default '\n\n')
            how to separate data from different TimeSeries objects
            (default is two blank lines)
        print_labels : bool {True, False}
            whether to print column labels on first line

        Returns
        -------
        str
        """
        experiment_path = self.parser.experiment.abspath
        if self.obs is None:
            raise ValueError('observable is not defined')

        printout = '# Data of plot {}\n#\n'.format(self.label)
        info = ('# observable:{}{}'.format(sep, self.obs.name) + '\n'
                '# experiment:{}{}'.format(sep, experiment_path) + '\n'
                '# filterset:{}{}'.format(sep, repr(self.parser.fset)) + '\n'
                '# condition set:{}{}'.format(sep, repr(self.cset)) + '\n')
        printout += info + '#\n'

        for index, (lin, ts) in enumerate(self._samples):
            if index == 0 and print_labels:
                local = ts.as_text(sep=sep, cell_sep=cell_sep,
                                   print_labels=print_labels)
            else:
                local = ts.as_text(sep=sep, cell_sep=cell_sep,
                                   print_labels=False)

            printout += local + timeseries_sep
        return printout.lstrip().rstrip()

    def save(self, user_bname=None, user_directory=None, extension='.pdf',
             add_obs=True, with_data_text=True, **kwargs):
        """Save figure (and data as text if specified).

        Parameters
        ----------
        user_directory : str (default None)
            directory under which plot and metadata are saved. When None,
            canonical path is sampleplots under experiment folder
        user_bname : str
            label of the file to be saved
        extension : str (default .pdf)
            will decide the format of output graphics
        add_obs : bool {True, False}
            append observable string representation to bname
        with_data_text : bool {True, False}
            whether to export a formatted text file with data
        kwargs : keyword arguments
            keyword arguments of self.data_as_text() method

        Note
        ----
        Call to save must succede to a make_plot call that defines obs.
        """
        if self.obs is None:
            raise ValueError('observable is not defined')
        exp = self.parser.experiment
        if user_directory is None:
            path = os.path.join(exp.abspath, 'sampleplots')
            if not os.path.exists(path):
                os.makedirs(path)
        else:
            path = user_directory
        if user_bname is None:
            bname = self.label
        else:
            bname = user_bname
        if add_obs:
            bname += '-' + self.obs.name
        figname = os.path.join(path, bname + '-plot' + extension)
        self.figpath = figname
        self.fig.savefig(figname, bbox_inches='tight')
        print('Figure saved as {}'.format(figname))
        if with_data_text:
            dataname = os.path.join(path, bname + '-data.txt')
            with open(dataname, 'w') as f:
                f.write(self.data_as_text(**kwargs))
        return


def plot_samples(samples, obs, parser=None, conditions=[],
                 report_condition='master',
                 units='',
                 yscale='linear',
                 show_markers=True,
                 marker='o',
                 show_lines=True,
                 linestyle='-',
                 join_cells=False,
                 end_points_emphasis=False,
                 color='C0',
                 change_cell_color=False,
                 change_lineage_color=False,
                 change_colony_color=False,
                 change_container_color=False,
                 alpha=1.,
                 superimpose='none',
                 limit_axes=20,
                 axe_xsize=6,
                 xrange_fractional_pad=.1,
                 axe_xrange=(None, None),  # auto
                 axe_ysize=1.6,
                 yrange_fractional_pad=.2,
                 axe_yrange=(None, None),  # auto
                 yrange_nticks=2,
                 report_cids=True,
                 report_cids_yposAxes=.8,
                 report_divisions=True,
                 show_legend=True,
                 data_statistics=False,
                 ref_mean=None, ref_var=None):
    """Plot small samples trajectories with various display options.

    Main options are the superimposing of timeseries on each subplot, the
    coloring of trajectories, and more.
    Ordering of samples matters: each sample is a couple (Lineage, TimeSeries)
    and lineages should be grouped by colony and containers. Coloring options
    only work well when samples are correctly ordered.

    Parameters
    ----------
    samples : list of (Lineage instance, TimeSeries instance)
        samples to be plotted
    obs : Observable instance
        observable to be plotted
    parser : Parser instance
        when data_statistics is used, needs parser.fset info
    conditions : list of :class:`FilterSet` instances
        set of conditions that are evaluated
    units : str (default '')
        units to add to title
    report_condition: str (default 'master')
        which condition to report on the plot: This must be either 'master',
        either the repr of one item of conditions, or its label attribute
    suppl_obs: list of :class:`Observables`
        supplementary observables that need to be computed so that each filter
        in conditions can be applied
    yscale : str {'linear', 'log'}
        scale of y-axis
    show_markers : bool {True, False}
        whether to plot markers for data (True is encouraged)
    marker : str
        type of markers
    show_lines : bool {True, False}
        whether to show lines. When show_markers is active, lines will get
        a low a low alpha value (0.3); when inactive, alpha parameters adjusts
        line transparency
    join_cells : bool {False, True}
        whether to connect parent last frame to daughter first frame with
        dashed line
    end_points_emphasis : bool {False, True}
        whether to add markers for first (circles) and last (squares) frames
        for each cell
    color : str {C0, C1, ..., C9}
        matplotlib range spec
    change_cell_color : bool {False, True}
        if active, cell color depends on cell generation index in colony
        priority over change_lineage_color
    change_lineage_color : bool {False, True}
        if active and change_cell_color inactive, each lineage gets a different
        color.
    change_colony_color : bool {False, True}
        when previous change_color options are False, set whether color changes
        from one colony to the next
    change_container_color : bool {False, True}
        when previous change_color options are False, set whether color changes
        from one container to the next
    alpha : float (default 1.)
        opacity of main artist: either markers if active, either lines
    superimpose : {'none', 'all', 'container', 'colony', int}
        where int is the number of lineages to superimpose per axe;
        when 'none' is equivalent to superimpose=1;
        when 'all', all lineages from colony are drawn on same axe; this is
        identical to setting superimpose=len(colony.idseqs);
        when 'related', all consecutive lineages that belongs to the same
        colony are plotted on the same axe
    limit_axes : int (default 50)
        limit the number of axes on a figure (will limit size as well, as each
        axe has fixed size)
    axe_xsize : float
        size of the x-axis (inches)
    xrange_fractional_pad : float (default .1)
        when xrange is fixed automatically by getting min, max values, set the
        padding to left, rightas a fraction of max-min range
    axe_xrange : couple of floats (default (None, None))
        user can force values for the left, right bounds here
    axe_ysize : float
        size if a single ax y-axis (inches)
    yrange_fractional_pad : float (default .2)
        when yrange is fixed automatically by getting min, max values, set the
        padding to bottom, top as a fraction of max-min range
    axe_yrange : couple of floats (default (None, None))
        user can force values for the bottom, top bounds here
    yrange_nticks : int (default 2)
        number of ticks on the y-axis
    report_cids : {False, True}
        whether to report for cell identifiers on top of axes
    report_cids_yposAxes : float (between 0 and 1, default .8)
        y position for reporting cell identifiers, in transAxes units
    report_divisions : {False, True}
        whether to draw vertical lines spanning multiple axes to visualize
        the tree structure arising from divisions
        this option is set to False automatically if superimpose is not 'none'
    show_legend : bool {True, False}
        whether to show legend
    data_statistics : bool {True, False} or str
        if True, try to find analysis folder to plot mean values from data;
        if str, tries to identify the appropriate condition from conditions
        to represent (can be 'master', repr(condition), or condition.label)
    ref_mean : float (default None)
        user can set the expected value here
    ref_var : float (default None)
        user can set the variance here
    """
    # check that report_condition is a valid condition
    condition_repr = 'master'
    if conditions:
        if report_condition in map(repr, conditions):
            condition_repr = report_condition  # default setting
            # find it
            for cdt in conditions:
                if repr(cdt) == report_condition:
                    condition_human_readable = str(cdt)
        # if one used the FilterSet.label
        else:
            for cdt in conditions:
                if report_condition == cdt.label:
                    condition_repr = repr(cdt)
                    condition_human_readable = report_condition
    if condition_repr == 'master':
        condition_human_readable = 'No condition'

    # counting the number of axes and makinf dictionary lineage index: ax index
    numbers = []
    nsuperimpose = 1  # default
    if superimpose == 'none':  # default setting to visualize tree structure
        n_axes = len(samples)
        if n_axes > limit_axes:
            msg = 'Too many lineages ({});'.format(len(samples))
            msg += ' reducing to first {}'.format(limit_axes)
            warnings.warn(msg)
            n_axes = limit_axes
        numbers = [1 for iax in range(n_axes)]
        nsuperimpose = 1
    elif superimpose == 'all':
        n_axes = 1
        nsuperimpose = len(samples)
        report_divisions = False
        numbers = [len(samples)]
    elif superimpose == 'colony':
        numbers = _count_colonies(samples)
        n_axes = len(numbers)
        if n_axes > limit_axes:
            msg = 'Too many lineages ({});'.format(len(samples))
            msg += ' reducing to first {}'.format(sum(numbers[:limit_axes]))
            warnings.warn(msg)
            n_axes = limit_axes
        nsuperimpose = max(numbers)
    elif superimpose == 'container':
        numbers = _count_containers(samples)
        n_axes = len(numbers)
        if n_axes > limit_axes:
            msg = 'Too many lineages ({});'.format(len(samples))
            msg += ' reducing to first {}'.format(sum(numbers[:limit_axes]))
            warnings.warn(msg)
            n_axes = limit_axes
        nsuperimpose = max(numbers)
    elif isinstance(superimpose, int):
        n_axes = len(samples) // superimpose
        if len(samples) % superimpose > 0:
            n_axes += 1
        if n_axes > limit_axes:
            new_lim = limit_axes * superimpose
            n_axes = limit_axes
            msg = 'Too many lineages ({});'.format(len(samples))
            msg += ' reducing to first {}'.format(new_lim)
            warnings.warn(msg)
        nsuperimpose = superimpose
        numbers = [nsuperimpose for iax in range(n_axes)]
    # make dictionary
#    print(numbers)
    index_to_iax = _make_dic(numbers)

    # if more than 1 timeseries per ax, discard reporting cids and divisions
    if nsuperimpose > 1:
        report_cids = False
        report_divisions = False

    # we will check at the end if each axe is not empty
    at_least_one_timeseries = np.array([False for iax in range(n_axes)])

    # setting figure and axes
    fig, axes = plt.subplots(n_axes, 1,
                             figsize=(axe_xsize, axe_ysize * n_axes))

    # compatibility 1 and >1 axes
    if n_axes == 1:
        axes = [axes, ]

    # storing min, max values for setting boundaries at the end
    lefts, rights, bottoms, tops = [], [], [], []

    # user-defined values
    if ref_mean is not None:
        for ax in axes:
            ax.axhline(ref_mean, ls='-.', color='C7', alpha=1.)
            if ref_var is not None:
                ax.axhline(ref_mean + np.sqrt(ref_var), ls=':', color='C7')
                ax.axhline(ref_mean - np.sqrt(ref_var), ls=':', color='C7')

    # plotting values from computed statistics
    data_stat_handles = []
    if data_statistics:
        statistics_condition_repr = 'master'  # default
        # data_statistics can be True, repr(condition), condition.label
        if isinstance(data_statistics, str):
            statistics_condition_repr = data_statistics
        if parser is None:
            raise IOError('Need parser to get fset info')
        else:
            try:

                res = add_data_statistics(axes, parser, obs, conditions,
                                          condition_repr=statistics_condition_repr)
                for item in res:
                    if item is not None:
                        data_stat_handles.append(item)
            except IOError:
                msg = ('Statistics have not been computed yet.')
                warnings.warn(msg)

    line2D_valid = None
    line2D_unvalid = None
    line2D_join = None

    container_lab = None
    colony_root = None
    iax = None

    pcolor = re.compile('C(\d)')  # pattern matching

    # looping through samples
    for index, (lineage, ts) in enumerate(samples):
        # let's pick the corresponding axe
        try:
            this_iax = index_to_iax[index]
        except KeyError:
            break  # we're out of registered axes
        # if we've reached the end of axes, let's get out of here
        if this_iax >= limit_axes:
            break

        ax = axes[this_iax]

        this_container = lineage.colony.container.label
        this_root = lineage.colony.root

        # priority for changing color: cell > lineage > colony > container
        if change_cell_color:
            # set color level depending on tree depth of first plotted cell
            cell = lineage.colony.get_node(lineage.idseq[0])
            depth = lineage.colony.depth(node=cell)
            color = 'C{}'.format(depth % 10)
        elif change_lineage_color:
            color = 'C{}'.format(index % 10)
        elif change_colony_color:
            # test when we change colony
            boo = (container_lab is not None and
                   (this_container != container_lab or
                    this_root != colony_root))
            if boo:
                # match previous
                m = pcolor.match(color)
                if m:
                    scindex, = m.groups()
                    cindex = int(scindex)
                    color = 'C{}'.format((cindex + 1) % 10)
        elif change_container_color:
            # test when we change colony
            boo = (container_lab is not None and
                   this_container != container_lab)
            if boo:
                # match previous
                m = pcolor.match(color)
                if m:
                    scindex, = m.groups()
                    cindex = int(scindex)
                    color = 'C{}'.format((cindex + 1) % 10)

        # write container label to each new panel
        show_colony_root = False  # needs activation
        if this_iax != iax:
            show_colony_root = (this_root != colony_root and
                                superimpose in ['none', 'colony', 1])
            msg = ''
            if show_colony_root:
                 msg += 'cont. {}, root {}'.format(this_container, this_root)
                 ax.text(0.01, 0.05, msg, color='C7', alpha=.8, transform=ax.transAxes)
        if len(ts.timeseries.clear) > 0:
            at_least_one_timeseries[this_iax] = True
            x = ts.timeseries.clear.x
            y = ts.timeseries.clear.y
            # only NaNs ? move along
            if np.isnan(x).all() or np.isnan(y).all():
                continue
            # calling add_timeseries
            ret = add_timeseries(ax, ts, condition_repr=condition_repr,
                                 end_points_emphasis=end_points_emphasis,
                                 show_markers=show_markers,
                                 marker=marker,
                                 show_lines=show_lines,
                                 linestyle=linestyle,
                                 join_cells=join_cells,
                                 color=color,
                                 alpha=alpha,
                                 change_cell_color=change_cell_color,
                                 use_last_color=True,
                                 report_cids=report_cids,
                                 report_cids_yposAxes=report_cids_yposAxes,
                                 )

            # getting min, max values; possible to get color: dat.get_color()
            limits, lines = ret
            left, right, bottom, top = limits
            lines_valid, lines_unvalid, lines_join = lines
            lefts.append(left)
            rights.append(right)
            bottoms.append(bottom)
            tops.append(top)

            if line2D_valid is None and lines_valid:
                line2D_valid = lines_valid[0]
            if line2D_unvalid is None and lines_unvalid:
                line2D_unvalid = lines_unvalid[0]
            if line2D_join is None and lines_join:
                line2D_join = lines_join[0]

        container_lab = this_container
        colony_root = this_root
        iax = this_iax

    # PLOT SETTINGS

    # found limits
    left = np.amin(lefts)
    right = np.amax(rights)
    bottom = np.amin(bottoms)
    top = np.amax(tops)

    hrange = right - left
    vrange = top - bottom
    hfrac = xrange_fractional_pad  # fraction of horizontal blank space to the leftmost, rightmost point
    vfrac = yrange_fractional_pad  # fraction of vertical blank space to the bottom, top point

    # user-defined limits
    if axe_xrange[0] is not None:
        if axe_xrange[0] > left:
            warnings.warn('Using non-inclusive user-defined left bound')
        left = axe_xrange[0]
    else:
        left = left - hfrac * hrange
    if axe_xrange[1] is not None:
        if axe_xrange[1] < right:
            warnings.warn('Using non-inclusive user-defined right bound')
        right = axe_xrange[1]
    else:
        right = right + hfrac * hrange
    if axe_yrange[0] is not None:
        if axe_xrange[0] > bottom:
            warnings.warn('Using non-inclusive user-defined bottom bound')
        bottom = axe_yrange[0]
    else:
        bottom = bottom - vfrac * vrange
    if axe_yrange[1] is not None:
        if axe_xrange[1] < top:
            warnings.warn('Using non-inclusive user-defined top bound')
        top = axe_yrange[1]
    else:
        top = top + vfrac * vrange

    # ticks
    # locator
    locator = _set_time_axis_ticks(axes[0], obs, bounds=(left, right))
    for ax in axes:
#        ax.xaxis.set_major_locator(ticker.MultipleLocator(60))
        ax.xaxis.set_major_locator(locator)

    # erase labels for intermediate layers
    for ax in axes[1:-1]:
        ax.set_xticklabels([])

    for iax, ax in enumerate(axes):
        ax.tick_params(axis='x', direction='in')
        ax.set_xlim(left=left, right=right)
        ax.set_ylim(bottom=bottom, top=top)
        if yscale == 'log':
            ax.set_yscale('log')
        if not at_least_one_timeseries[iax]:
            ax.text(0.4, 0.4, "NO DATA", transform=ax.transAxes)

    if len(axes) > 1:
        for ax in axes[:-1]:
            ax.spines['bottom'].set_visible(False)
            ax.tick_params(axis='x', colors='C7')
        for ax in axes[1:]:
            ax.spines['top'].set_color('C7')

        axes[0].tick_params(axis='x', direction='in', bottom='on',
                            labelbottom='off')

    yfmt = ticker.ScalarFormatter()
    yfmt.set_powerlimits((-1, 3))

    ax = axes[0]
    loc = ticker.MaxNLocator(nbins=yrange_nticks)
    ax.yaxis.set_major_locator(loc)
    # call loc to set up
    _ = loc()
    ax.yaxis.set_major_formatter(yfmt)
    yfmt.set_locs(ax.get_yticks())
    yticklocs = yfmt.locs
#    print(yticklocs)
    yticklabels = [yfmt.pprint_val(item) for item in yticklocs]

    # use first ax locations and formetter
    for ax in axes[1:]:
        ax.yaxis.set_major_locator(ticker.FixedLocator(yticklocs))
        ax.yaxis.set_ticklabels(yticklabels)

    axes[-1].set_xlabel('Time (mins)', x=.95, horizontalalignment='right',
                        fontsize='medium')

    # now that limits have been set : report for divisions
    if report_divisions:
        line2D_join = _report_divisions(samples, axes, line2D_join, index_to_iax, limit_axes)

    # add legend
    if show_legend:
        ax = axes[-1]
        clab = condition_human_readable
        handles = []
        labels = []
        if line2D_valid is not None:
            valid = mlines.Line2D([], [])
            valid.update_from(line2D_valid)
            lab_from_data = valid.get_label()
            color = valid.get_color()
            output_lab = 'tracks'
            if report_cids:
                output_lab += ' (e.g. cell {})'.format(lab_from_data)
            if report_condition != 'master':
                output_lab += ', {}: YES'.format(clab)
            valid.set_label(output_lab)
            valid.set_color(color)
            handles.append(valid)
            labels.append(valid.get_label())
        if line2D_unvalid is not None:
            unvalid = mlines.Line2D([], [])
            unvalid.update_from(line2D_unvalid)
            lab_from_data = unvalid.get_label()
            color = unvalid.get_color()
            output_lab = 'tracks'
            if report_cids:
                output_lab += ' (e.g. cell {})'.format(lab_from_data)
            if report_condition != 'master':
                output_lab += ', {}: NO'.format(clab)
            unvalid.set_label(output_lab)
            unvalid.set_color(color)
            handles.append(unvalid)
            labels.append(unvalid.get_label())
        if report_divisions and line2D_join is not None:
            join = mlines.Line2D([], [])
            join.update_from(line2D_join)
            color = join.get_color()
            lab = 'connections at division'
            join.set_label(lab)
            join.set_color(color)
            handles.append(join)
            labels.append(join.get_label())
        data_stat_labs = [item.get_label() for item in data_stat_handles]
        handles += data_stat_handles
        labels += data_stat_labs
        ax.legend(handles=handles, labels=labels,
                  loc='upper left',
                  bbox_to_anchor=(0, -.5/axe_ysize),
                  borderaxespad=0.)

    # add title
    titling = r'{}'.format(obs.as_latex_string)
    if units:
        titling += ' ({})'.format(units)
    axes[0].text(0.5, 1 + .2/axe_ysize, titling,
                size='large',
                horizontalalignment='center',
                verticalalignment='bottom',
                transform=axes[0].transAxes)

    fig.subplots_adjust(hspace=0)

    return fig


def _report_divisions(samples, axes, line2D_join, index_to_iax, limit_axes):
    """Helper function to draw vertical lines at cell divisions

    As this function use data, and axe inverted transforms, it needs to be
    called after all axe limits have been set

    Parameters
    ----------
    samples : list of (Lineage instance, corresponding TimeSeries
        samples to be plotted
    axes : list of Axes instances
        axes used to plot successive samples
    line2D_join : Line2D instance
        sets the styling when not None
    index_to_iax : dict
        maps sample index to axe index
    limit_axes : int
        maximum number of axes

    Returns
    -------
    line2D_join : Line2D instance
        used to get styling in main function
    """
    color = 'C7'  # uniform color for tree visualization
    vjoin = None
    for index, (lineage, ts) in enumerate(samples):
        # let's pick the corresponding axe
        try:
            this_iax = index_to_iax[index]
        except KeyError:
            break  # we're out of registered axes
        # if we've reached the end of axes, let's get out of here
        if this_iax >= limit_axes:
            break

        if len(ts.timeseries.clear) > 0:

            x = ts.timeseries.clear.x
            y = ts.timeseries.clear.y

            # visualisation of cell divisions, works only when superimpose='none'
            logging.debug('reporting cell division')
            # getting parameters
            if line2D_join is not None:
                jlw = line2D_join.get_linewidth()
                jls = line2D_join.get_linestyle()
                jalpha = line2D_join.get_alpha()
            # default
            else:
                jlw = 1.
                jls = ':'
                jalpha = .8
            # styling options dict
            styling = {'lw': jlw, 'ls': jls, 'color': color, 'alpha': jalpha}
            # get first cell
            if len(lineage.idseq) > 0:
                cid = lineage.idseq[0]
                cell = lineage.colony.get_node(cid)
                if cell.birth_time is not None:
                    pid = cell.bpointer
                    pvals = None  # must be set to Coordinates instances if found
                    # get index corresponding to pid AND last parent value
                    for pi, (plin, pts) in enumerate(samples[:index+1]):
                        # check that colonies match
                        if plin.colony == lineage.colony:
                            # parse previous lineage idseq to find pid
                            for pindex, plinid in enumerate(plin.idseq):
                                # get parent values
                                if pid == plinid:
                                    parent_iax = index_to_iax[pi]
                                    logging.debug('cid: {} (in axe {}), pid: {} (in axe {})'.format(cid, this_iax,pid, parent_iax))
                                    pslice = pts.slices[pindex]
#                                    x_pvals = pts.timeseries.x[pts.slices[pindex]]
#                                    y_pvals = pts.timeseries.y[pts.slices[pindex]]
#                                    pvals = Coordinates(x_pvals, y_pvals, x_name=pts.timeseries.x_name, y_name=pts.timeseries.y_name)
                                    pvals = pts.timeseries[pslice].clear
                                    break
                            if pid in plin.idseq:
                                break
                    # it works only when pid is in a previous lineage
                    if pi < index:
                        logging.debug('current lineage index {}, parent lineage index {}'.format(index, pi))
                        i_top = index_to_iax[pi]
                        # get transData and inverse transAxe
                        # last valid parent y value in transAxes coords
                        if len(pvals) > 0:
                            dispx, dispy = axes[i_top].transData.transform((pvals.x[-1], pvals.y[-1]))
                            inv = axes[i_top].transAxes.inverted()
                            _, ymax = inv.transform((dispx, dispy))
                        # when no valid values, point to center y
                        else:
                            ymax = .5  # in transAxes coordinates
                        # logging.debug('in parent cell axe, ymax={}'.format(ymax))
                        vjoin = axes[i_top].axvline(cell.birth_time, ymax=ymax, **styling)
                        i_bottom = index_to_iax[index]
                        # first valid cell y value in transAxes coords
                        if len(x) > 0:
                            dispx, dispy = axes[i_bottom].transData.transform((x[0], y[0]))
                            inv = axes[i_bottom].transAxes.inverted()
                            _, ymin = inv.transform((dispx, dispy))
                        # when no valid values, point to center y
                        else:
                            ymin = .5  # in transAxes coords
                        # logging.debug('in current cell axe, ymin={}'.format(ymin))
                        axes[i_bottom].axvline(cell.birth_time, ymin=ymin, **styling)

                        for k in range(i_top + 1, i_bottom):
                            axes[k].axvline(cell.birth_time, **styling)
    if line2D_join is None and vjoin is not None:
        line2D_join = vjoin
    return line2D_join


def _count_containers(samples):
    label = None
    lin_per_container = []
    count_lin = 0
    for lin, ts in samples:
        this_label = lin.colony.container.label
        if label is not None and this_label != label:
            lin_per_container.append(count_lin)
            count_lin = 0
        label = this_label
        count_lin += 1
    if count_lin != 0:
        lin_per_container.append(count_lin)
    return lin_per_container


def _count_colonies(samples):
    label = None
    croot = None
    count_lin = 0
    lin_per_colony = []
    for lin, ts in samples:
        this_label = lin.colony.container.label
        this_root = lin.colony.root
        # changing container
        if label is not None and this_label != label:
            lin_per_colony.append(count_lin)
            count_lin = 0
        elif croot is not None and this_root != croot:
            lin_per_colony.append(count_lin)
            count_lin = 0
        croot = this_root
        label = this_label
        count_lin += 1
    if count_lin != 0:
        lin_per_colony.append(count_lin)
    return lin_per_colony


def _make_dic(numbers):
    """Make dictionary lineage index to axe index when there are numbers[k]
    lineages for axe k.

    Parameters
    ----------
    numbers : list of int
        size is size of plot

    Returns
    -------
    dictionary: int to int
        key: lineage index
        value: axe index
    """
    dic = {}
    cumul = np.cumsum(numbers)
    for k in range(len(cumul)):
        if k == 0:
            sl = range(cumul[k])
        else:
            sl = range(cumul[k-1], cumul[k])
        for item in sl:
            dic[item] = k
    return dic
