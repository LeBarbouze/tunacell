#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module defines filters for Cell instances
"""
from __future__ import print_function

import copy
import numpy as np

from tunacell.filters.main import FilterGeneral, bounded, included, FilterAND
from tunacell.base.datatools import multiplicative_increments
from tunacell.base.observable import Observable


class FilterCell(FilterGeneral):
    "General class for filtering cell objects (reader.Cell instances)"

    _type = 'CELL'


class FilterCellAny(FilterCell):
    "Class that does not filter anything."

    def __init__(self):
        self.label = 'Always True'  # short description for human readers
        return

    def func(self, cell):
        return True


class FilterData(FilterCell):
    """Default filter test only if cell exists and cell.data non empty."""

    def __init__(self):
        self.label = 'Cell Has Data'
        return

    def func(self, cell):
        boo = False
        if cell is not None:
            boo = cell.data is not None and len(cell.data) > 0
        return boo


class FilterCellIDparity(FilterCell):
    """Test whether identifier is odd or even"""

    def __init__(self, parity='even'):
        self.parity = parity
        self.label = 'Cell identifier is {}'.format(parity)
        return

    def func(self, cell):
        # test if even
        try:
            even = int(cell.identifier) % 2 == 0
            if self.parity == 'even':
                return even
            elif self.parity == 'odd':
                return not even
            else:
                raise ValueError("Parity must be 'even' or 'odd'")
        except ValueError as ve:
            print(ve)
            return False


class FilterCellIDbound(FilterCell):
    """Test class"""

    def __init__(self, lower_bound=None, upper_bound=None):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        self.label = '{} <= cellID < {}'.format(lower_bound, upper_bound)
        return

    def func(self, cell):
        return bounded(int(cell.identifier),
                       self.lower_bound, self.upper_bound)


class FilterHasParent(FilterCell):
    """Test whether a cell has an identified parent cell"""

    def __init__(self):
        self.label = 'Cell Has Parent'
        return

    def func(self, cell):
        boo = False
        if cell.parent:
            boo = True
        return boo


class FilterDaughters(FilterCell):
    "Test whether a given cell as at least one daughter cell"

    def __init__(self, daughter_min=1, daughter_max=2):
        label = 'Number of daughter cell(s): '
        label += '{0} <= n_daughters <= {1}'.format(daughter_min, daughter_max)
        self.label = label
        self.lower_bound = daughter_min
        self.upper_bound = daughter_max + 1  # lower <= x < upper
        return

    def func(self, cell):
        return bounded(len(cell.childs),
                       lower_bound=self.lower_bound,
                       upper_bound=self.upper_bound)


class FilterCompleteCycle(FilterCell):
    "Test whether a cell has a given parent and at least one daughter."

    def __init__(self, daughter_min=1):
        label = 'Cell cycle complete'
        label += ' (with at least {} daughter cell(s)'.format(daughter_min)
        self.daughter_min = daughter_min
        self.label = label
        return

    def func(self, cell):
        filt_parent = FilterHasParent()
        filt_daughter = FilterDaughters(daughter_min=self.daughter_min)
        return filt_parent(cell) and filt_daughter(cell)


class FilterCycleFrames(FilterCell):
    """Check whether cell has got a minimal number of datapoints."""

    def __init__(self, lower_bound=None, upper_bound=None):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        label = 'Number of registered frames:'
        label += '{0} <= n_frames <= {1}'.format(self.lower_bound,
                                                 self.upper_bound)
        self.label = label
        return

    def func(self, cell):
        # check whether data exists
        boo = False
        filtData = FilterData()
        if filtData.func(cell):
            boo = bounded(len(cell.data),
                          lower_bound=self.lower_bound,
                          upper_bound=self.upper_bound
                          )
        return boo


class FilterCycleSpanIncluded(FilterCell):
    """Check that cell cycle time interval is within valid bounds."""

    def __init__(self, lower_bound=None, upper_bound=None):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        label = '{} <= Initial frame and Final frame < {}'.format(lower_bound,
                                                                  upper_bound)
        self.label = label
        return

    def func(self, cell):
        boo = False
        filtData = FilterData()
        if filtData(cell):
            boo = included(cell.data['time'],
                           lower_bound=self.lower_bound,
                           upper_bound=self.upper_bound)
        return boo


class FilterTimeInCycle(FilterCell):
    """Check that tref is within cell birth and division time"""

    def __init__(self, tref=0.):
        self.tref = tref
        label = 'birth/first time <= {} < division/last time'.format(tref)
        self.label = label
        return

    def func(self, cell):
        boo = False
        filtData = FilterData()
        if filtData(cell):
            if cell.birth_time is not None:
                lower = cell.birth_time
            else:
                lower = cell.data['time'][0]
            if cell.division_time is not None:
                upper = cell.division_time
            else:
                upper = cell.data['time'][-1]
            boo = lower <= self.tref < upper
        return boo


class FilterObservableBound(FilterCell):
    """Check that a given observable is bounded.

    Parameters
    ----------
    obs : Observable instance
        observable that will be tested for bounds
        works only for continuous observable (mode='dynamics')
    tref : float (default None)
        Time of reference at which to test dynamics observable value
    lower_bound : float (default None)
    upper_bound : float (default None)
    """

    def __init__(self, obs=Observable(name='undefined'), tref=None,
                 lower_bound=None, upper_bound=None):
        self.obs_to_test = obs  # observable to be tested
        self._obs = [obs, ]  # hidden to be computed at for filtering purpose
        self.tref = tref
        # code below is commented because func is able to deal with arrays
#        if obs.mode == 'dynamics' and tref is None:
#            msg = 'For dynamics mode, this filter needs a valid tref (float)'
#            raise ValueError(msg)
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        label = '{} <= {}'.format(lower_bound, obs.name)
        if tref is not None:
            label += ' (t={})'.format(tref)
        label += ' < {}'.format(upper_bound)
        self.label = label
        return

    def func(self, cell):
        import collections
        boo = False
        if self.tref is not None:
            filt = FilterAND(FilterData(),
                             FilterTimeInCycle(tref=self.tref))
        else:
            filt = FilterData()
        label = self.obs_to_test.label
        if filt(cell):
            # retrieve data
            array = cell._sdata[label]  # two cases: array, or single value
            raw_time = cell.data['time']
            if len(raw_time) > 1:
                dt = np.amin(raw_time[1:] - raw_time[:-1])
            else:
                dt = cell.container.period
            if array is None:
                return False
            if isinstance(array, collections.Iterable):
                if self.tref is None:
                    # data may be one value (for cycle observables), or array
                    boo = bounded(array[label], self.lower_bound, self.upper_bound)
                else:
                    # find data closest to tref (-> round to closest time)
                    # for now return closest time to tref
                    index = np.argmin(np.abs(array['time'] - self.tref))
                    # check that it's really close:
                    if np.abs(array['time'][index] - self.tref) < dt:
                        value = array[label][index]
                        boo = bounded(value, self.lower_bound, self.upper_bound)
            # otherwise it's a number
            else:
                boo = bounded(array, self.lower_bound, self.upper_bound)
        return boo


# useless ?
class FilterLengthIncrement(FilterCell):
    "Check increments are bounded."

    def __init__(self, lower_bound=None, upper_bound=None):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        label = 'Length increments between two successive frames: '
        label += '{0} <= delta_length <= {1}'.format(self.lower_bound,
                                                     self.upper_bound)
        return

    def func(self, cell):
        boo = False
        filtData = FilterData()
        if filtData(cell):
            ell = np.array(cell.data['length'])
            incr = multiplicative_increments(ell)
            lower = bounded(np.amin(incr), lower_bound=self.lower_bound)
            upper = bounded(np.amax(incr), upper_bound=self.upper_bound)
            boo = lower and upper
        return boo


class FilterSymmetricDivision(FilterCell):
    """Check that cell division is (roughly) symmetric.

    Parameters
    ----------
    raw : str
        column label of raw observable to test for symmetric division
        (usually one of 'length', 'area'). This quantity will be approximated
    """

    def __init__(self, raw='area', lower_bound=0.4, upper_bound=0.6):
        self.raw = raw
        # Observable to be computed: raw at birth, raw at division
        # hidden _obs because not part of parameters, but should be computed
        self._obs = [Observable(raw=raw, scale='log', mode='birth', timing='b'),
                     Observable(raw=raw, scale='log', mode='division',
                                timing='d')]
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        label = 'Symmetric division filter:'
        ratio_str = '(daughter cell {})/(mother cell {})'.format(raw, raw)
        label += ' {} <= {} <= {}'.format(self.lower_bound,
                                          ratio_str,
                                          self.upper_bound)
#        label += 'OR (in case mother cell data is missing) '
#        label += '{0} <= (daughter cell area)/(sister cell area) <= {1}\
#                  '.format(self.lower_bound/self.upper_bound,
#                           self.upper_bound/self.lower_bound)
        self.label = label
        return

    def func(self, cell):
        boo = False
        filtData = FilterData()
        if cell.parent is None:
            # birth is not reported, impossible to test, cannot exclude from data
            boo = True
        else:
            if filtData(cell):
                csize = cell._sdata[self._obs[0].label]
                if filtData(cell.parent):
                    psize = cell.parent._sdata[self._obs[1].label]
                    boo = bounded(csize/psize,
                                  lower_bound=self.lower_bound,
                                  upper_bound=self.upper_bound
                                  )
                else:
                    # parent exists, but without data.
                    # this is a weird scenario, that should not exist
                    # TODO: warn user?
                    # but we can check with sibling
                    sibs = copy.deepcopy(cell.parent.childs)
                    for item in sibs:
                        if item.identifier == cell.identifier:
                            sibs.remove(item)
                    if sibs:
                        if len(sibs) > 1:
                            from ..base.cell import CellChildsError
                            raise CellChildsError('>2 daughters')
                        sib = sibs[0]  # there should be only one cell
                        if sib.data is not None:
                            sibsize = sib._sdata[self._obs[0].label()]
                            boo = bounded(csize/sibsize,
                                          lower_bound=self.lower_bound,
                                          upper_bound=self.upper_bound
                                          )
                        else:
                            boo = True  # sibling cell: no data, accept this cell
                    else:
                        boo = True  # no sibling, accept this cell
        return boo
