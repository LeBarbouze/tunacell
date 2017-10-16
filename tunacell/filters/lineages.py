#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Classes to filter lineages.
"""
from __future__ import print_function


import numpy as np

from tunacell.filters.main import FilterGeneral, bounded, intersect
from tunacell.filters.cells import FilterData

from io import StringIO


class FilterLineage(FilterGeneral):
    "TODO: define lineage object (as tree?)"

    _type = 'LINEAGE'


class FilterLineageAny(FilterLineage):
    "No selection"

    def __init__(self):
        self.label = 'Always True'
        return

    def func(self, lineage):
        return True


class FilterLineageData(FilterLineage):
    "Select lineages with non-empty data"

    def __init__(self):
        self.label = 'Lineage with non-empty data'
        return

    def func(self, lineage):
        boo = False
        if lineage is not None:
            if lineage.data is not None:
                if len(lineage.data) > 0:
                    boo = True
        return boo


class FilterLineageLength(FilterLineage):
    "Select lineages with bounded number of cells"

    def __init__(self, lower_bound=1, upper_bound=None):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.label = '{} <= number of cells <= {}'.format(lower_bound,
                                                          upper_bound)
        return

    def func(self, lineage):
        if bounded(len(lineage.idseq),
                   lower_bound=self.lower_bound,
                   upper_bound=self.upper_bound):
            return True
        else:
            return False


class FilterLineageTimeIntersect(FilterLineage):
    """Select lineages which data time interval has non-empty intersection
    with given bounds.
    """
    def __init__(self, lower_bound=None, upper_bound=None):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        label = 'At least one time frame between {} and {}'.format(lower_bound,
                                                                   upper_bound)
        self.label = label
        return

    def func(self, lineage):
        boo = False
        filtData = FilterLineageData()
        if filtData(lineage):
            boo = intersect(lineage.data['time'],
                            lower_bound=self.lower_bound,
                            upper_bound=self.upper_bound)
        return boo


class FilterLineageTimeBound(FilterLineage):
    """Select lineages which data time interval is bounded by parameters"""

    def __init__(self, lower_bound=None, upper_bound=None):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        label = '{} < lineage time range < {}'.format(lower_bound, upper_bound)
        self.label = label
        return

    def func(self, lineage):
        boo = False
        filtData = FilterLineageData()
        if filtData(lineage):
            boo = bounded(lineage.data['time'],
                          lower_bound=self.lower_bound,
                          upper_bound=self.upper_bound)
        return boo


class FilterLineageTimeLength(FilterLineage):
    """Select Lineages which data time interval is bounded by given params."""

    def __init__(self, lower_bound=0., upper_bound=None):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.label = '{} <= time span <= {}'.format(lower_bound, upper_bound)
        return

    def func(self, lineage):
        times = lineage.data['time']
        tmin = np.amin(times)
        tmax = np.amax(times)
        if bounded(tmax - tmin,
                   lower_bound=self.lower_bound,
                   upper_bound=self.upper_bound):
            return True
        else:
            return False


# %% Next Filters are used for conditional analysis on lineages

class FilterLineageWithCellProperty(FilterLineage):
    """Initialize instance

    Parameters
    ----------
    cell_filter : Filter instance acting on Cells
        must be callable with a Cell instance as parameter
    extend_ancestry: boolean
        whether to look only to cells forming lineage sequence (False),
        or to look to all ancestor cells til colony root (True)
    """

    def __init__(self, cell_filter=FilterData(), extend_ancestry=True):
        self.cell_filter = cell_filter
        self._obs.extend(cell_filter._obs)  # add hidden observables
        self.extend_ancestry = extend_ancestry
        label = 'A cell in lineage '
        if extend_ancestry:
            label += '(up to root) '
        label += 'with:\n'
        f = StringIO(u'{}'.format(cell_filter))  # StringIO needs unicode
        for line in f.readlines():
            label += '    {}'.format(str(line))  # decode to appropriate str
        label += '\n'
        self.label = label.rstrip()
        return

    def func(self, lineage):
        boo = False
        if len(lineage.idseq) > 0:
            leaf = lineage.idseq[-1]
            for cid in lineage.colony.rsearch(leaf):
                # if we don't extend to root, stop when we're out of lineage
                if not self.extend_ancestry:
                    if cid not in lineage.idseq:
                        break

                cell = lineage.colony.get_node(cid)
                if self.cell_filter(cell):
                    boo = True
                    break  # no need to go further
        return boo
