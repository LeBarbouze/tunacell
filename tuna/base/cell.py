#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module defines how cells are stored as tuna's objects
"""
from __future__ import print_function

import numpy as np
import warnings
from copy import deepcopy

import treelib as tlib


from tuna.datatools import (local_rate, extrapolate_endpoints,
                            derivative, logderivative, ExtrapolationError)


class CellError(Exception):
    pass


class CellChildsError(CellError):
    pass


class CellParentError(CellError):
    pass


class CellDivisionError(CellError):
    pass


class Cell(tlib.Node):
    """General class to handle cell data structure.

    Inherits from treelib.Node class to facilitate tree building.

    Parameters
    ----------
    identifier : str
        cell identifier
    container : :class:`Container` instance
        container to which cell belongs

    Attributes
    ----------
    container : :class:`Container` instance
        container to chich cell belongs
    childs : list of :class:`Cell` instances
        daughter cells of current cell
    parent : :class:`Cell` instance
        mother cell of current cell
    birth_time : float (default None)
        time of cell birth (needs to be computed)
    division_time : float (default None)
        time of cell division (needs to be computed)

    Methods
    -------
    set_division_events()
        computes birth/division times when possible
    build(obs)
        builds timeseries, uses one of the following methods depending on obs
    build_timelapse(obs)
        builds and stores timeseries associated to obs, in 'dynamics' mode
    build_cyclized(obs)
        builds and stores cell-cycle value associated to obs, not in 'dynamics'
        mode
    """

    def __init__(self, identifier=None, container=None):

        tlib.Node.__init__(self, identifier=identifier)

        self._childs = []
        self._parent = None
        self._birth_time = None
        self._division_time = None
        self._sdata = {}  # dictionary to contain computed data
        self.container = container  # point to Container instance
        # cells are built from a specific container instance
        # container can be a given field of view, a channel, a microcolony, ...

        return

    # We add few definitions to be able to chain between Cell instances
    @property
    def childs(self):
        "Get list of child instances."
        return self._childs

    @childs.setter
    def childs(self, value):
        if value is None:
            self._childs = []
        elif isinstance(value, list):
            for item in value:
                self.childs = item
        elif isinstance(value, Cell):
            self._childs.append(value)
        else:
            raise CellChildsError

    @property
    def parent(self):
        "Get parent instance."
        return self._parent

    @parent.setter
    def parent(self, pcell):
        if pcell is None:
            self._parent = None
        elif isinstance(pcell, Cell):
            self._parent = pcell
        else:
            raise CellParentError

    @property
    def birth_time(self):
        "Get cell cycle start time. See below for Setter."
        return self._birth_time

    @birth_time.setter
    def birth_time(self, value):
        "Set cell cycle start time. See above for Getter."
        self._birth_time = value

    @property
    def division_time(self):
        "Get cell cycle end time. See below for Setter."
        return self._division_time

    @division_time.setter
    def division_time(self, value):
        "Set cell cycle end time. See above for Getter."
        if self.birth_time is not None:
            if value < self.birth_time:
                raise CellDivisionError
        self._division_time = value

    def set_division_event(self):
        "method to call when parent is identified"
        previous_frame = None
        if (self.parent is not None) and (self.parent.data is not None):
            previous_frame = self.parent.data['time'][-1]

        first_frame = None
        if self.data is not None:
            first_frame = self.data['time'][0]

        if previous_frame is not None and first_frame is not None:
            div_time = (previous_frame + first_frame)/2.  # halfway
            self.birth_time = div_time
            self.parent.division_time = div_time

        return

    def __repr__(self):
        cid = str(self.identifier)
        if self.parent:
            pid = str(self.parent.identifier)
        else:
            pid = '-'
        if self.childs:
            ch = ','.join([c.identifier for c in self.childs])
        else:
            ch = '-'
        return cid+';p:'+pid+';ch:'+ch

    def info(self):
        dic = {}
        dic['a. Identifier'] = '{}'.format(self.identifier)
        pid = 'None'
        if self.parent:
            pid = '{}'.format(self.parent.identifier)
        dic['b. Parent id'] = pid
        chids = 'None'
        if self.childs:
            chids = ', '.join(['{}'.format(ch.identifier)
                               for ch in self.childs])
        dic['c. Childs'] = chids
        dic['d. Birth time'] = '{}'.format(self.birth_time)
        dic['e. Division time'] = '{}'.format(self.division_time)
        if self.data is not None:
            dic['f. N_frames'] = '{}'.format(len(self.data))
        return dic

    def build(self, obs):
        """Builds timeseries"""
        if obs.mode == 'dynamics':
            return self.build_timelapse(obs)
        else:
            return self.compute_cyclized(obs)

    def build_timelapse(self, obs):
        """Builds timeseries corresponding to observable of mode 'dynamics'.

        Parameters
        ----------
        obs : Observable instance
            mode must be 'dynamics'

        Returns
        -------
        Numpy structured array with 3 columns 'time', <obs.label>, 'cellID'

        Notes
        -----
        Some observables carry the 'local_fit' option True. In this case,
        local fits over shifting time-windows are performed. If one would keep
        only a given cell's data, then the constraints on shifting time-window
        would let some 'empty' times, at which no evaluation can be performed.
        This is solved by getting data from the cell's parent cell's data. This
        operation computes time-window fiited data in the cell's parent cycle.
        Two precautions must then be taken:
            1. a given cell's data must be used only once for evaluating parent
               cell's data,
            2. when data has been used from one daughter cell, concatenate
               the current cell's evaluated data to it.

        .. warning::
           For some computations, the time interval between consecutive
           acquisitions is needed. If it's defined in the container or the
           experiment metadata, this parameter will be imported; otherwise if
           there are at least 2 consecutive values, it will be inferred from
           data (at the risk of making mistakes if there are too many missing
           values)
        """
        label = str(obs.label())
        yaxis = obs.raw
        # if empty, return empty array of appropriate type
        if len(self.data) == 0:  # there is no data, but it has some dtype
            arr = np.array([], dtype=[('time', 'f8'),
                                      (label, 'f8'),
                                      ('cellID', self.data.dtype['cellID'])])
            return arr  # empty array is returned with appropriate dtype

        dt = self.container.period
        if dt is None:
            # automatically finds dt
            if len(self.data) > 1:
                arr = self.data['time']
                time_increments = arr[1:] - arr[:-1]
                dt = np.round(np.amin(np.abs(time_increments)), decimals=2)

        # define which function to apply to cell to retrieve individual
        # timeseries
        # case: differentiation
        time, array = [], []
        if obs.differentiate:

            # case : no local fit, use finite differences
            if not obs.local_fit:
                if obs.scale == 'linear':
                    if len(self.data) > 1:
                        out = derivative(self.data[['time', yaxis]])
                        time, array = out
                elif obs.scale == 'log':
                    if len(self.data) > 1:
                        out = logderivative(self.data[['time', yaxis]])
                        time, array = out
                idarray = self.data['cellID'][:len(time)]  # resize
            # case : local estimates using local_rates
            else:
                # dismiss the adjusted values
                fit, adjusted = local_rate(self, yaxis=yaxis,
                                           yscale=obs.scale,
                                           time_window=obs.time_window,
                                           dt=dt,
                                           join_points=obs.join_points)
                time = fit['time']
                array = fit['rate_' + yaxis]
                idarray = fit['cellID']

        # case: no differentiation
        else:
            # case : no local fitting, retrieve data
            if not obs.local_fit:
                time = self.data['time']
                array = self.data[yaxis]
                idarray = self.data['cellID']
            # case : local fits using local_rate
            else:
                fit, adjusted = local_rate(self, yaxis=yaxis,
                                           yscale=obs.scale,
                                           time_window=obs.time_window,
                                           dt=dt,
                                           join_points=obs.join_points)
                time = fit['time']
                array = fit['fit_' + yaxis]
                idarray = fit['cellID']

        # when local fitting is used, the result may span 2 cells:
        # current cell, and its parent
        out = np.zeros(len(time), dtype=[('time', 'f8'),
                                         (label, 'f8'),
                                         ('cellID', idarray.dtype)])
        out['time'] = time[:]
        out[label] = array[:]
        out['cellID'] = idarray[:]

        # update cell values
        if idarray.dtype.kind in ['i', 'u']:
            boo = idarray == int(self.identifier)  # text idtype is integer
        else:
            boo = idarray == self.identifier  # text type is set as strings

        # check if non-empty and non-intersecting
        if label in self._sdata.keys():
            to_concat = out[boo]
            existing = self._sdata[label]
            if _disjoint_time_sets(to_concat['time'], existing['time']):
                self._sdata[label] = np.concatenate((to_concat, existing))
                # this should be well sorted since operations on parents
                # would provide data towards end of cell-cycle
        else:
            self._sdata[label] = out[boo]

        # update parent cell values if existing and not updated yet by sibling
        if self.parent is not None:
            pid = self.parent.identifier
            if idarray.dtype.kind in ['i', 'u']:
                boo = idarray == int(pid)
            else:
                boo = idarray == pid
            new = out[boo]
            # not updated yet
            if label not in self.parent._sdata.keys():
                self.parent._sdata[label] = new
            # check that present data is time disjoint
            elif _disjoint_time_sets(new['time'],
                                     self.parent._sdata[label]['time']):
                arr = self.parent._sdata[label]
                self.parent._sdata[label] = np.concatenate((arr, new))

#        # update parent cell values if they were not updated by sibling
#        if self.bpointer is not None:
#            # check that sibling has not been computed
#            fill_parent = False
#            # TODO ERASE TRY BLOCK WHEN UNBUGGED
##            try:
##                if len(self.parent.childs) > 1:
##                    pass
##            except AttributeError as ae:
##                print(self.identifier)
##                raise ae
#            ###
#            if len(self.parent.childs) > 1:
#                for ch in self.parent.childs:
#                    if ch.identifier != self.identifier:
#                        break
#                if label not in ch._sdata.keys():
#                    fill_parent = True
#            else:
#                fill_parent = True
#            if fill_parent:
#                boo = idarray == int(self.parent.identifier)
#                new = out[boo]
#                if label not in self.parent._sdata.keys():
#                    self.parent._sdata[label] = new
#                else:
#                    # collect previously computed values, concatenate new
#                    arr = self.parent._sdata[label]
#                    self.parent._sdata[label] = np.concatenate((arr, new))
        return out

    def compute_cyclized(self, obs):
        """Computes observable when mode is different from 'dynamics'.

        Parameters
        ----------
        obs : Observable instance
            mode must be different from 'dynamics'

        Returns
        -------
        float corresponding to desired observable

        Raises
        ------
        ValueError
            when Observable mode is 'dynamics'
        """
        scale = obs.scale
        npts = obs.join_points
        label = obs.label()
        if obs.mode == 'dynamics':
            raise ValueError('Called build_cyclized for dynamics mode')
        # associate continous observable and build corresponding ._sdata
        cobs = deepcopy(obs)
        cobs.mode = 'dynamics'
        cobs.timing = 't'
        clabel = cobs.label()
        # discard result as it can mix cell, and parent cell data
        _ = self.build_timelapse(cobs)
        # now we compute cell cycle observable using created _sdata: only cell
        time = self._sdata[clabel]['time']
        array = self._sdata[clabel][clabel]
        # get value
        try:
            if obs.mode == 'birth':
                value = extrapolate_endpoints(self,
                                              zip(time, array),
                                              scale=scale,
                                              end_point='birth',
                                              join_points=npts)
            elif obs.mode == 'division':
                value = extrapolate_endpoints(self,
                                              zip(time, array),
                                              scale=scale,
                                              end_point='division',
                                              join_points=npts)
            elif 'net-increase' in obs.mode:
                dval = extrapolate_endpoints(self,
                                             zip(time, array),
                                             scale=scale,
                                             end_point='division',
                                             join_points=npts)
                bval = extrapolate_endpoints(self,
                                             zip(time, array),
                                             scale=scale,
                                             end_point='birth',
                                             join_points=npts)
                if obs.mode == 'net-increase-additive':
                    value = dval - bval
                elif obs.mode == 'net-increase-multiplicative':
                    value = dval/bval
            elif obs.mode == 'average':
                value = np.nanmean(array)
            elif obs.mode == 'rate':
                if obs.scale == 'log':
                    array = np.log(array)
                value, intercept = np.polyfit(time, array, 1)
        except ExtrapolationError as err:
            msg = '{}'.format(err)
            warnings.warn(msg)
            value = np.nan  # missing information
        self._sdata[label] = value
        return value


def _disjoint_time_sets(ts1, ts2):
    if len(ts1) == 0 or len(ts2) == 0:
        return True
    min1, min2 = map(np.nanmin, [ts1, ts2])
    max1, max2 = map(np.nanmax, [ts1, ts2])
    return max1 < min2 or max2 < min1
