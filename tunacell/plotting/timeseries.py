#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module provides few functions to ease plotting TimeSeries objects
"""
from __future__ import print_function

import numpy as np
import re

import matplotlib as mpl
import matplotlib.transforms as transforms

from tunacell.stats.api import load_univariate, UnivariateIOError

def add_data_statistics(axes, parser, obs, conditions,
                        condition_repr='master'):
    """Add dynamics of single observable statistics.

    Parameters
    ----------
    axes : list of :class:`matplotlib.axes.Axes`
    parser : :class:`Parser` instance
    obs : :class:`Observable` instance
    conditions : list of :class:`FilterSet` instances
    condition_repr : str (default 'master')
        use computed statistics for this condition;
        must be either 'master', either a repr of one of conditions.selections

    Returns
    -------
    line_mean, fill_std
    line_mean : matplotlib.Lines2D
    fill_std matplotlib.collections.PolyCollection
    """
    try:
        single = load_univariate(parser.experiment, obs, cset=conditions)
    except UnivariateIOError:
        return None, None
    # condition_repr can be either a repr(condition), condition.label, or 'master'
    use_repr = 'master'
    human_readable = 'all'
    if condition_repr != 'master':
        for cdt in conditions:
            if condition_repr == repr(cdt):
                use_repr = condition_repr
                human_readable = cdt.label
                break
            elif condition_repr == cdt.label:
                use_repr = repr(cdt)
                human_readable = condition_repr
                break
    item = single[use_repr]
    tt = item.time
    mm = item.average
    std = item.std

    for ax in axes:
        line_mean, = ax.plot(tt, mm, color='C7', alpha=0.8,
                             label='average ({} samples)'.format(human_readable))
        fill_std = ax.fill_between(tt, mm - std, mm + std,
                                   facecolor='C7', alpha=0.2,
                                   label='+/-1 standard deviation')
    return line_mean, fill_std


def add_timeseries(ax, ts, condition_repr='master',
                   show_markers=True,
                   marker='o',
                   end_points_emphasis=False,
                   show_lines=True,
                   linestyle='-',
                   join_cells=False,
                   color=None,
                   alpha=1.,
                   change_cell_color=False,
                   use_last_color=True,
                   report_cids=False,
                   report_cids_yposAxes=.9,
                   ):
    """Plot timeseries object in ax.

    Parameters
    ----------
    ax : :class:`matplotlib.axes.Axes` instance
        the axis onto which plotting is performed
    ts: :class:`tunacell.base.timeseries.Timeseries` instance
        TimeSeries objects contain data
    condition_repr : str (default 'master')
        repr of :class:`FilterSet` instance used to condition data, or 'master'
    show_markers: bool {True, False}
        whether to report data as data points (*i.e.* with markers)
    marker : str
        marker style to use for data-points
    end_points_emphasis : bool {False, True}
        whether to plot extra, larger markers for first and last frames
    show_lines : bool {True, False}
        whether to report data as lines. When active, check whether
        show_markers is active: display line with reduced transparency,
        otherwise use alpha paramater below
    linestyle : str
        linestyle to use
    join_cells : bool {False, True}
        whether to plot a dashed line connected parent to daughter cell
    color : str (default None)
        color to use to plot the first cell in the TimeSeries
    alpha : float (default 1.)
        transparency with which points/lines are plotted.
        Note that when both are plotted, lines have a reduced transparency.
    change_cell_color : bool {False, True}
        whether to change color from one cell to the next
    use_last_color : bool {True, False}
        whether to use last line color plotted on the current Axe
    report_cids : bool {True, False}
        whether to report cell identifiers above datapoints
    report_cids_yposAxes : float (default .9)
        yaxis position to print cell identifiers, in Axe coordinates (stay
        between 0 and 1 to be in the vertical space defined by current Axe)

    Returns
    -------
    (left, right, bottom, top), (lines_valid, lines_unvalid, joins)

    (left, right, bottom, top): floats
       boundaries values from data plotted (useful to determine ax.ylim - xlim)
    lines_valid : list of Line2D instances
        only main Line2D are reported (when both markers and lines are plotted
        only the Line2D object associated to markers is reported here)
        Those ones are the timeseries that fulfils condition_repr
    lines_unvalid: list of Line2D instances
        idem, but for timeseries that does not fulfil condition_repr
    joins : list of Line2D instances
        lines connecting cells when join_cells is True
    """
    # size settings defined by matplotlib.rcParams 
    markersize = mpl.rcParams['lines.markersize']
    markeredgewidth = mpl.rcParams['lines.markeredgewidth']
    linewidth = mpl.rcParams['lines.linewidth']

    # PLOTTING PARAMETERS
    alpha_connecting = .8
    if not show_markers and not show_lines:
        show_lines = True
    # color to be used
    ucolor = None  # last used color (undefined at first call)
    ncolor = None  # color for the next cell (when change_cell_color is True)
    if color is not None:
        ucolor = color
    # pattern matching
    pcolor = re.compile('C(\d)')

    # find boundaries on unconditioned timeseries
    x = ts.timeseries.x
    y = ts.timeseries.y
    left, bottom = map(np.nanmin, [x, y])
    right, top = map(np.nanmax, [x, y])

    if len(x) > 1:
        xstep = np.nanmin(x[1:]-x[:-1])
    else:
        xstep = 5.

    inits, fins = [], []
    prej = None
    # define heterogeneous transform
    # the x coords of this transformation are data, and the
    # y coord are axes
    trans = transforms.blended_transform_factory(ax.transData, ax.transAxes)

    # we plot the timeseries cell by cell
    lines_valid = []
    lines_unvalid = []
    joins = []

    dat = None  # Line2D for data points
    connecting = None  # Line2D for connecting lines
    for si, sl in enumerate(ts.slices):
        # get cell identifier
        cid = ts.ids[si]
        # check that sl is a proper slice
        if not isinstance(sl, slice):
            # it means no data for this cell
            if report_cids:
                # check that there were previous data and print, otherwise do nothing
                if len(fins) > 0:
                    xposData = fins[-1][0]
                    ax.text(xposData + xstep, report_cids_yposAxes,
                            '{}'.format(cid),
                            transform=trans, color='k', alpha=.5)
            continue

        # check that data is valid under condition_repr
        if not ts.selections[condition_repr][si]:
            appears = False
        else:
            appears = True

        # get data for corresponding slice (i.e. corresponding cell)
        local = ts.timeseries.x[sl]
        if len(local) > 0:
            # plot each cell data
            # data points
            xdata = ts.timeseries[sl].clear.x
            ydata = ts.timeseries[sl].clear.y
            # is there valid data?
            # if empty, move to next slice
            if len(xdata) == 0:
                # update color is using change_cell_color
                if change_cell_color and color is not None:
                    m = pcolor.match(color)
                    if m:
                        sindex, = m.groups()
                        index = int(sindex)
                        ncolor = 'C{}'.format((index+1) % 10)
                continue
            # COLOR MANAGEMENT
            # if we keep same color, keep same color as previous cell
            if not change_cell_color:
                color = ucolor
            elif ncolor is not None:
                color = ncolor
            if color is not None:
                m = pcolor.match(color)
                if m:
                    sindex, = m.groups()
                    index = int(sindex)
                    ncolor = 'C{}'.format((index+1) % 10)
            dat, = ax.plot(xdata, ydata,
                           ls='None',
                           marker=marker,
                           markersize=markersize,
                           markeredgewidth=markeredgewidth,
                           color=color,
                           alpha=alpha,
                           label='{}'.format(cid)
                           )
            ucolor = dat.get_color()  # use identical color
            if not appears:
                dat.set_markerfacecolor('None')
            if not show_markers:
                dat.set_visible(False)
            # first/last points
            if end_points_emphasis:
                x, y = xdata[0], ydata[0]
                firstframe, = ax.plot(x, y, marker='s',
                                      markersize=markersize + 2.,
                                      markeredgewidth=markeredgewidth,
                                      markerfacecolor='None',
                                      markeredgecolor=ucolor,
                                      label='{}: first frame'.format(cid)
                                      )
                x, y = xdata[-1], ydata[-1]
                lastframe, = ax.plot(x, y, marker='o',
                                     markersize=markersize + 2.,
                                     markeredgewidth=markeredgewidth,
                                     markerfacecolor='None',
                                     markeredgecolor=ucolor,
                                     label='{}: last frame'.format(cid)
                                     )
            # lines
            connecting, = ax.plot(xdata, ydata,
                                  marker=None,
                                  ls=linestyle,
                                  lw=linewidth,
                                  color=ucolor,
                                  label='connecting {}'.format(cid)
                                  )
            if not appears:
                connecting.set_linestyle('--')
            if not show_lines:
                connecting.set_visible(False)
            elif show_markers:
                connecting.set_alpha(alpha_connecting)
            else:
                connecting.set_alpha(alpha)  # no data points: use full alpha

            if show_markers:
                if appears:
                    lines_valid.append(dat)
                else:
                    lines_unvalid.append(dat)
            else:
                if appears:
                    lines_valid.append(connecting)
                else:
                    lines_unvalid.append(connecting)
            if report_cids:
                xposData = np.percentile(xdata, 35.)
                ax.text(xposData, report_cids_yposAxes, '{}'.format(cid),
                        transform=trans, color=color, alpha=.8)
            # append first/last frame (registered if cell completed its cycle)
            if sl.start is not None:
                inits.append((xdata[0], ydata[0]))
                if join_cells and prej is not None and sl.start == prej:  # ensure consecutive
                    xjoins, yjoins = zip(*(fins[-1:] + inits[-1:]))
                    join, = ax.plot(xjoins, yjoins,
                                    ls=':',
                                    lw=linewidth,
                                    alpha=alpha_connecting,
                                    color=ucolor)
                    joins.append(join)
#                    if not join_cells:
#                        join.set_visible(False)
            if sl.stop is not None:
                fins.append((xdata[-1], ydata[-1]))
        # there might be no data for current slice, but still the cell exists
        else:
            if report_cids:
                if len(fins) > 0:
                    xposData = fins[-1][0]
                    ax.text(xposData + xstep, report_cids_yposAxes,
                            '{}'.format(cid),
                            transform=trans, color='k', alpha=.5,
#                            fontsize=fontsize
                            )
        prej = sl.stop

    return (left, right, bottom, top), (lines_valid, lines_unvalid, joins)
