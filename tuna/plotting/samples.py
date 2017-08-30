#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module defines how samp
"""
from __future__ import print_function

import os
import warnings
import re
import datetime
import inspect

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import collections

from tuna.base.colony import Colony
from tuna.base.lineage import Lineage
from tuna.plotting.timeseries import add_timeseries, add_data_statistics

try:
    from StringIO import StringIO  # python2
except ImportError:
    from io import StringIO  # python3


class SamplePlot(object):
    """General class to store Colony plots.

    Parameters
    ----------
    obs : :class:`Observable` instance
        observable to be plotted
    samples : {(list of) :class:`Colony` instance(s),
               (list of) :class:`Lineage` instance(s)}
        timeseries from these samples will be plotted
    parser : :class:`Parser` instance
    conditions : list of :class:`FilterSet` instances

    """

    def __init__(self, obs, samples, parser=None, conditions=[]):
        self.obs = obs
        self.parser = parser
        self.cset = conditions

        # building sample list
        self._samples = []
        if isinstance(samples, collections.Iterable):
            for sample in samples:
                if isinstance(sample, Colony):
                    for lin in sample.iter_lineages():
                        self._add_atomic_sample(lin)
                elif isinstance(sample, Lineage):
                    self._add_atomic_sample(sample)
        elif isinstance(samples, Colony):
            for lin in samples.iter_lineages():
                self._add_atomic_sample(lin)
        elif isinstance(samples, Lineage):
            self._add_atomic_sample(samples)

        today = datetime.datetime.today()
        self.label = 's{}-{}'.format(today.strftime('%Y%m%d%H%M%S'), str(obs))

        self.fig = None  # Figure instance

        return

    def make_plot(self, report_condition='master',
                  yscale='linear',
                  show_markers=True,
                  marker='o',
                  markersize=6.,
                  markeredgewidth=.8,
                  show_lines=True,
                  linestyle='-',
                  linewidth=2.,
                  join_cells=True,
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
                  axe_ysize=1.6,
                  report_cids=True,
                  report_cids_yposAxes=.8,
                  report_divisions=True,
                  show_legend=True,
                  data_statistics=False,
                  ref_mean=None, ref_var=None):
        """Produce Figure and save plotted data.

        Look at plot_samples doctring for parameters (too long to repeat here)
        """
        # make keyword arguments
        plt.ioff()
        kwargs = {}
        frame = inspect.currentframe()
        res = inspect.getargvalues(frame)
        for key in res.args:
            if key == 'self':
                continue
            kwargs[key] = res.locals[key]
        # make figure by calling plot_samples
        fig = plot_samples(self._samples, self.obs,
                           conditions=self.cset,
                           parser=self.parser,
                           **kwargs)
        self.fig = fig
        return

    def _add_atomic_sample(self, lineage):
        ts = lineage.get_timeseries(self.obs, self.cset)
        self._samples.append((lineage, ts))
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

        printout = '# Data of plot {}\n#\n'.format(self.label)
        info = ('# experiment:{}{}'.format(sep, experiment_path) + '\n'
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
        """
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
                bname += '-' + str(self.obs)
        figname = os.path.join(path, bname + '-plot' + extension)
        self.fig.savefig(figname, bbox_inches='tight')
        if with_data_text:
            dataname = os.path.join(path, bname + '-data.txt')
            with open(dataname, 'w') as f:
                f.write(self.data_as_text(**kwargs))
        return


def plot_samples(samples, obs, parser=None, conditions=[],
                 report_condition='master',
                 yscale='linear',
                 show_markers=True,
                 marker='o',
                 markersize=6.,
                 markeredgewidth=.8,
                 show_lines=True,
                 linestyle='-',
                 linewidth=2.,
                 join_cells=True,
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
                 axe_ysize=1.6,
                 fontsize=10.,
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
    samples : list of (Lineage instance, corresponding TimeSeries
        samples to be plotted
    obs : Observable instance
        observable to be plotted
    parser : Parser instance
        when data_statistics is used, needs parser.fset info
    conditions : list of :class:`FilterSet` instances
        set of conditions that are evaluated
    report_condition: str (default 'master')
        which condition to report on the plot: This must be either 'master',
        either a repr of one element of conditions
    suppl_obs: list of :class:`Observables`
        supplementary observables that need to be computed so that each filter
        in conditions can be applied
    yscale : str {'linear', 'log'}
        scale of y-axis
    show_markers : bool {True, False}
        whether to plot markers for data (True is encouraged)
    marker : str
        type of markers
    markersize : float
        size of markers
    markeredgewidth : float
        width of marker edges; used for cell note verifying condition_label
    show_lines : bool {True, False}
        whether to show lines. When show_markers is active, lines will get
        a low a low alpha value (0.3); when inactive, alpha parameters adjusts
        line transparency
    join_cells : bool {True, False}
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
        size of the x-axis
    axe_ysize : float
        size if a single ax y-axis
    fontsize : float
        basis of fontsizes,in points
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
    data_statistics : bool {True, False}
        if True, try to find analysis folder to plot mean values from data
    ref_mean : float (default None)
        user can set the expected value here
    ref_var : float (default None)
        user can set the variance here
    """
    # get fontsize
    default_fs = mpl.rcParams['font.size']
    mpl.rc('font', size=fontsize)
    # check that report_condition is a valid condition
    condition_label = 'master'
    if conditions:
        if report_condition in map(repr, conditions):
            condition_label = report_condition  # default setting
            # find it
            for cdt in conditions:
                if repr(cdt) == report_condition:
                    condition_human_readable = str(cdt)
    if condition_label == 'master':
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
        if parser is None:
            raise IOError('Need parser to get fset info')
        else:
            try:
                res = add_data_statistics(axes, parser, obs, conditions,
                                          condition_label=condition_label)
                for item in res:
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
        show_container_label = False  # needs activation
        show_colony_root = False  # needs activation
        if this_iax != iax:
            show_container_label = (this_container != container_lab and
                                    (superimpose in ['none',
                                                     'colony',
                                                     'container',
                                                     1]))
            if show_container_label:
                msg = 'Container {}'.format(this_container)
                ax.text(0.01, 0.05, msg, color='C7', transform=ax.transAxes)
            show_colony_root = (this_root != colony_root and
                                superimpose in ['none', 'colony', 1])
            if show_colony_root:
                msg = 'Root {}'.format(this_root)
                ax.text(0.01, 0.9, msg, color=color, transform=ax.transAxes)

        if len(ts.timeseries.clear_x) > 0:
            at_least_one_timeseries[this_iax] = True
            x = ts.timeseries.clear_x
            y = ts.timeseries.clear_y
            # only NaNs ? move along
            if np.isnan(x).all() or np.isnan(y).all():
                continue
            # calling add_timeseries
            ret = add_timeseries(ax, ts, condition_label=condition_label,
                                 end_points_emphasis=end_points_emphasis,
                                 show_markers=show_markers,
                                 marker=marker,
                                 markersize=markersize,
                                 markeredgewidth=markeredgewidth,
                                 show_lines=show_lines,
                                 linestyle=linestyle,
                                 linewidth=linewidth,
                                 join_cells=join_cells,
                                 color=color,
                                 alpha=alpha,
                                 change_cell_color=change_cell_color,
                                 use_last_color=True,
                                 report_cids=report_cids,
                                 report_cids_yposAxes=report_cids_yposAxes,
                                 fontsize='medium')

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

            # visualisation of cell divisions, works only when superimpose='none'
            if report_divisions:
                # getting parameters
                if line2D_join is not None:
                    jlw = line2D_join.get_linewidth()
                    jls = line2D_join.get_linestyle()
                    jalpha = line2D_join.get_alpha()
                # default
                else:
                    jlw = 1.
                    jls = '--'
                    jalpha = .3
                # get first cell
                if len(lineage.idseq) > 0:
                    cid = lineage.idseq[0]
#                    print(cid)
                    cell = lineage.colony.get_node(cid)
                    if cell.birth_time is not None:
                        # add cell_birth as left limit
                        lefts.append(cell.birth_time)
                        pid = cell.bpointer
#                        print(pid)
                        # get index corresponding to pid
                        for pi, (plin, pts) in enumerate(samples[:index+1]):
                            if pid in plin.idseq and plin.colony == lineage.colony:
                                break
                        if pi < index:
                            i_top = index_to_iax[pi]
                            axes[i_top].axvline(cell.birth_time,
                                                ymax=.5,
                                                lw=jlw,
                                                ls=jls,
                                                color=color,
                                                alpha=jalpha)
                            i_bottom = index_to_iax[index]
                            axes[i_bottom].axvline(cell.birth_time, ymin=.5,
                                                   lw=jlw,
                                                   ls=jls,
                                                   color=color,
                                                   alpha=jalpha)
                            for k in range(i_top + 1, i_bottom):
                                axes[k].axvline(cell.birth_time,
                                                lw=jlw,
                                                ls=jls,
                                                color=color,
                                                alpha=jalpha)
        container_lab = this_container
        colony_root = this_root

    # PLOT SETTINGS
    left = np.amin(lefts)
    right = np.amax(rights)
    bottom = np.amin(bottoms)
    top = np.amax(tops)

    hrange = right - left
    vrange = top - bottom

    # ticks
    for ax in axes[1:-1]:
        ax.set_xticklabels([])

    for iax, ax in enumerate(axes):
        ax.tick_params(axis='x', direction='in')
        ax.set_xlim(left=left-hrange/10., right=right+hrange/10.)
        ax.set_ylim(bottom=bottom-vrange/10., top=top+vrange/10.)
        if yscale == 'log':
            ax.set_yscale('log')
        if not at_least_one_timeseries[iax]:
            ax.text(0.4, 0.4, "NO DATA", transform=ax.transAxes)

    axes[0].tick_params(axis='x', direction='out', top='on',
                        labeltop='on')
    axes[0].tick_params(axis='x', direction='in', bottom='on',
                        labelbottom='off')

    axes[-1].set_xlabel('Time (mins)', x=.95, horizontalalignment='right',
                        fontsize='large')
    if n_axes > 1:
        axes[0].xaxis.set_label_position('top')
        axes[0].set_xlabel('Time (mins)', x=.95, horizontalalignment='right',
                           fontsize='large')

    # add legend
    if show_legend:
        ax = axes[-1]
        color = 'C7'
        clab = condition_human_readable
        handles = []
        labels = []
        if line2D_valid is not None:
            valid = mlines.Line2D([], [])
            valid.update_from(line2D_valid)
            lab_from_data = valid.get_label()
            output_lab = ('data (e.g. cell {})'.format(lab_from_data) + ', '
                          'valid: {}'.format(clab))
            valid.set_label(output_lab)
            valid.set_color(color)
            handles.append(valid)
            labels.append(valid.get_label())
        if line2D_unvalid is not None:
            unvalid = mlines.Line2D([], [])
            unvalid.update_from(line2D_unvalid)
            lab_from_data = unvalid.get_label()
            output_lab = ('data (e.g. cell {})'.format(lab_from_data) + ', '
                          'not valid: {}'.format(clab))
            unvalid.set_label(output_lab)
            unvalid.set_color(color)
            handles.append(unvalid)
            labels.append(unvalid.get_label())
        if report_divisions and line2D_join is not None:
            join = mlines.Line2D([], [])
            join.update_from(line2D_join)
            lab = 'connections at division'
            join.set_label(lab)
            join.set_color(color)
            handles.append(join)
            labels.append(join.get_label())
        data_stat_labs = [item.get_label() for item in data_stat_handles]
        handles += data_stat_handles
        labels += data_stat_labs
        ax.legend(handles=handles, labels=labels,
                  loc=2,
                  bbox_to_anchor=(0, -.4),
                  borderaxespad=0.)

    # add title
    ax = axes[0]
    ax.text(0.5, 1.3, r'{}'.format(obs.as_latex_string),
            size='x-large',
            horizontalalignment='center',
            verticalalignment='bottom',
            transform=ax.transAxes)

    fig.subplots_adjust(hspace=0)
    # restore default font size
    mpl.rc('font', size=default_fs)
    return fig


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
