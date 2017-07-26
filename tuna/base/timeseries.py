#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module will define useful objects for conditional analysis
"""
import numpy as np
import pandas as pd


# define an object to handle heterogeneous types of time series
class TimeSeries(object):
    """Object that decorates the data with other useful attributes.

    Parameters
    ----------
    label: str
       label to be given to yaxis values (e.g. 'dot_width')
       this label may be used for subsequent output labelling
    ts : sequence of couples (time, value) (or numpy array)
    ids : sequence of cell identifiers from which data was collected
    index_cycles : sequence of couples (index_first, index_last)
       that delimit data corresponding to cell id
    select_ids : sequences of True/False values corresponding whether or
       not to include data from cell id in timeseries
    """

    def __init__(self, label, ts=[], ids=[], index_cycles=[], slices=None,
                 time_bounds=[],
                 select_ids={}):
        self.label = label  # label is the string label to be given to yaxis
#        self._groups = groups
        self._timeseries = ts
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
        label : str (default 'master')
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
        List of couples (time, value) for valid cells
        """
        selection = self.selections[condition_label]
        toconcat = []
        for index, cid in enumerate(self.ids):
            if selection[index] and self.slices[index] is not None:
                if sharp_tleft is not None:
                    if self.time_bounds[index][0] < sharp_tleft:
                        continue
                if sharp_tright is not None:
                    if self.time_bounds[index][1] > sharp_tright:
                        continue
                toconcat.append(self.timeseries[self.slices[index]])
        if len(toconcat) > 0:
            out = np.concatenate(toconcat)
        else:
            out = np.array([], dtype=[('time', 'f8'), (self.label, 'f8')])
        if len(out) > 0:
            if 'time' in out.dtype.names and self.label in out.dtype.names:
                out = out[['time', self.label]]
        return out

    @property
    def timeseries(self):
        return self._timeseries

    @timeseries.setter
    def timeseries(self, ts):
        self._timeseries = ts

    def __getitem__(self, key):
        return self.timeseries[key]

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
        labels = None
        if isinstance(self.timeseries, np.ndarray):
            labels = self.timeseries.dtype.names
        printout = ''
        if print_labels and labels is not None:
            printout += '\t'.join(labels) + '\n'
        printout += '\n'
        for sl in self.slices:
            chunk = ''
            local = self.timeseries[sl]
            for line in local:
                chunk += '{}'.format(sep).join(['{}'.format(item) for item in line]) + '\n'
            printout += chunk
            printout += cell_sep
        return printout.lstrip().rstrip()  # remove empty lines at beginning/end

    def to_dataframe(self, start_index=0):
        dic = {}
        # initialize dict of empty lists
        dic['time'] = []
        dic[self.label] = []
        dic['id'] = []
        for label in self.selections.keys():
            dic[label] = []
        for index, cid in enumerate(self.ids):
            if self.slices[index] is not None:
                ts = self.timeseries[self.slices[index]]
                times = ts['time']
                dic['time'].extend(times)
                dic[self.label].extend(ts[self.label])
                dic['id'].extend([cid for _ in times])
                for label in self.selections.keys():
                    boo = self.selections[label][index]
                    dic[label].extend([boo for _ in times])
        df = pd.DataFrame(dic, index=range(start_index, start_index + len(dic['time'])))
        return df
