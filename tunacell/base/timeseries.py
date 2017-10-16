#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module will define useful objects for conditional analysis
"""
import collections
import numpy as np
import pandas as pd

from tunacell.base.datatools import Coordinates


# define an object to handle heterogeneous types of time series
class TimeSeries(object):
    """Object that decorates the data with other useful attributes.

    Parameters
    ----------
    ts : :class:`Coordinates` instance, or 2d structured ndarray
        better to use Coordinates, so that names can be carried along
    ids : sequence of cell identifiers from which data was collected
    index_cycles : sequence of couples (index_first, index_last)
       that delimit data corresponding to cell id, must be same length as ids
    slices : sequence of slice objects
        each item can be used to slice the entire table
    time_bounds : sequence of couples of floats
        for each cell, first element is the lower bound of cell cycle, the
        second element is the upper bound of cell cycle, must be same length
        as ids
    select_ids : sequences of True/False values corresponding whether or
       not to include data from cell id in timeseries, must be same length as
       ids
    """

    def __init__(self, ts=[], ids=[], index_cycles=[], slices=None,
                 time_bounds=[], select_ids={}, container_label=None,
                 experiment_label=None):
        # ts is a Coordinates instance
        self.container_label = container_label
        self.experiment_label = experiment_label
        if isinstance(ts, Coordinates):
            self._timeseries = ts
        # ts is a numpy array (structured if possible)
        elif isinstance(ts, np.ndarray):
            # convert structured arrays to 2d ndarrays
            if ts.dtype.names is not None:
                _arr = ts.view((float, len(ts.dtype.names)))
                _x_name, _y_name = ts.dtype.names[:2]  # take only first 2 cols
            else:
                _arr = ts
                _x_name, _y_name = 'x', 'y'
            _x = _arr[:, 0]
            _y = _arr[:, 1]
            self._timeseries = Coordinates(_x, _y,
                                           x_name=_x_name, y_name=_y_name)
        # ... list of couples
        elif isinstance(ts, collections.Iterable):
            _ts = list(ts)
            _x, _y = map(np.array, zip(*_ts))
            self._timeseries = Coordinates(_x, _y)
        self.time_bounds = time_bounds
        self.slices = []
        if index_cycles:  # array indices corresponding to (first, last) frame for each cell
            self.index_cycles = index_cycles
            slices = []
            for item in index_cycles:
                if item is None:
                    slices.append(None)
                    # indices are reported as a single None
                    # when no data is reported for a given cell
                else:
                    i, j = item
                    if j is not None:
                        slices.append(slice(i, j+1))
                    else:
                        slices.append(slice(i, None))
            self.slices = slices
        elif slices is not None:
            self.slices = slices
            index_cycles = []
            for item in slices:
                if item is None:
                    index_cycles.append(None)
                else:
                    if item.stop is not None:
                        index_cycles.append((item.start, item.stop - 1))
                    else:
                        index_cycles.append((item.start, None))
            self.index_cycles = index_cycles
        self.ids = ids
        if len(select_ids.keys()) > 0:  # master is already defined
            self.selections = select_ids
        else:  # nothing is defined, we define master here
            self.selections = {'master': [True for _ in self.ids]}
        return

    def use_condition(self, condition_label='master',
                      sharp_tleft=None, sharp_tright=None):
        """Get conditioned timeseries.

        Parameter
        ---------
        condition_label : str (default 'master')
            must be a key of dictionary self.selections, and corresponds to
            the repr of a given :class:`FilterSet` instance.
        sharp_left : float (default None)
            sharp lower bound for cell cycle timing. USE ONLY FOR CELL CYCLE
            OBSERVABLES
        sharp_right : float (default None)
            sharp upper bound for cell cycle timing. USE ONLY FOR CELL CYCLE
            OBSERVABLES

        Returns
        -------
        Coordinates instance made of valid (x, y) points
        """
        selection = self.selections[condition_label]
        xs, ys = [], []
        for index, cid in enumerate(self.ids):
            if selection[index] and self.slices[index] is not None:
                if sharp_tleft is not None:
                    if self.time_bounds[index][0] < sharp_tleft:
                        continue
                if sharp_tright is not None:
                    if self.time_bounds[index][1] > sharp_tright:
                        continue
                xs.append(self.timeseries.x[self.slices[index]])
                ys.append(self.timeseries.y[self.slices[index]])
        if len(xs) > 0:
            _x = np.concatenate(xs)
            _y = np.concatenate(ys)
        else:
            _x = []
            _y = []
        out = Coordinates(_x, _y, x_name=self.timeseries.x_name,
                          y_name=self.timeseries.y_name)
        return out

    @property
    def timeseries(self):
        return self._timeseries
#
#    @timeseries.setter
#    def timeseries(self, ts):
#        self._timeseries = ts

#    def __getitem__(self, key):
#        return self.timeseries[key]

    def __repr__(self):
        return repr(self.timeseries)

    def as_text(self, sep='\t', cell_sep='\n', print_labels=False):
        """Export TimeSeries as text arrays

        Parameters
        ----------
        sep : str (default '\t')
            how to separate columns
        cell_sep : str (default '\n')
            how to separate cells (default: one blank line)
        print_labels : bool {False, True}
            first line is labels, followed by empty line
        """
        printout = ''
        labels = [self.timeseries.x_name,
                  self.timeseries.y_name,
                  'cellID',
                  'containerID',
                  'experimentID']
        if print_labels and labels is not None:
            printout += '\t'.join(labels) + '\n'
        printout += '\n'
        for index, sl in enumerate(self.slices):
            chunk = ''
            x = self.timeseries.x[sl]
            y = self.timeseries.y[sl]
            ids = len(x) * [self.ids[index]]
            container_id = len(x) * [self.container_label, ]
            exp_id = len(x) * [self.experiment_label, ]
            for line in zip(x, y, ids, container_id, exp_id):
                chunk += '{}'.format(sep).join(['{}'.format(item) for item in line]) + '\n'
            printout += chunk
            printout += cell_sep
        return printout.lstrip().rstrip()  # remove empty lines at beginning/end

    def to_dataframe(self, start_index=0, sharp_tleft=None, sharp_tright=None):
        dic = {}
        dic[self.timeseries.x_name] = []  # self.timeseries.x
        dic[self.timeseries.y_name] = []  # self.timeseries.y
        dic['cellID'] = []
        dic['containerID'] = []
        dic['experimentID'] = []
        for key in self.selections.keys():
            if key == 'master':
                continue
            dic[key] = []
        size = 0
        # add cell ID, container ID, experiment ID, and TRUE/FALSE for each cdt
        for index, sl in enumerate(self.slices):
            # collect only if within bounds
            if sharp_tleft is not None:
                if self.time_bounds[index][0] < sharp_tleft:
                    continue
            if sharp_tright is not None:
                if self.time_bounds[index][1] > sharp_tright:
                    continue
            _x = self.timeseries.x[sl]
            _y = self.timeseries.y[sl]
            dic[self.timeseries.x_name].extend(_x)
            dic[self.timeseries.y_name].extend(_y)
            dic['cellID'].extend(len(_x) * [self.ids[index], ])
            dic['containerID'].extend(len(_x) * [self.container_label, ])
            dic['experimentID'].extend(len(_x) * [self.experiment_label, ])
            # True/False for each
            for key, values in self.selections.items():
                # master: all True, useless to printout
                if key == 'master':
                    continue
                val = values[index]
                dic[key].extend(len(_x) * [val, ])
            size += len(_x)
        df = pd.DataFrame(dic, index=range(start_index, start_index + size))
        return df
