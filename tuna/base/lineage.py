#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module defines how to handle lineages.
"""
from __future__ import print_function

import numpy as np
import random
import collections

from tuna.observable import Observable

from tuna.base.timeseries import TimeSeries


class LineageError(Exception):
    pass


class LineageTimeseriesError(Exception):
    pass


class NoAncestry(Exception):
    pass


class Lineage(object):
    """General class to handle lineage data.

    Suitable for various conditioning methods.

    Parameters
    ----------
    tree : :class:`Colony` instance
    identifier_sequence : sequence of cell identifiers
        this is the (ordered) sequence of cells that defines the lineage

    Attributes
    ----------
    colony : :class:`Colony` instance
        used to derive lineage
    idseq : list of identifiers

    Methods
    -------
    get_generations(tref=None)
        Returns generation indices of cell sequence.
    get_boolean_tests(cset=None)
        Performs boolean test over cells in the sequence
    get_timeseries(obs, cset=None, suppl_obs=[], dt=5.)
        Returns :class:`Timeseries` instance of :class:`Observable`
        :param:`obs` in current lineage

    """

    def __init__(self, tree, identifier_sequence):
        self.colony = tree
        self.idseq = identifier_sequence
        self.division_timings = get_division_timing(self.idseq, tree)
        self._data = get_lineage_data(self.idseq, tree)
        # self.slices = self.split()
        return

    @property
    def data(self):
        return self._data

    def iter_cells(self, shuffle=False):
        """Iterator over cells in current lineage."""
        idseq = self.idseq
        if shuffle:
            idseq = random.shuffle(idseq)
        for cid in idseq:
            yield self.colony.get_node(cid)
        return

    def get_generations(self, tref=None):
        """Get generation number

        Parameter
        ---------
        tref : float (default None)
            optional to fix reference for generation count

        Returns
        -------
        list of int
            generation label within colony when tref=None, generation index
            since generation 0 at t=tref.

        Raises
        ------
        NoAncestry
            when tref is provided and no ancestry crosses tref.
        """
        gens = []
        genref = 0
        for cid in self.idseq:
            gens.append(self.colony.level(cid))
        gens = np.array(gens, dtype='i4')
        if tref is not None:
            from tuna.filters.cells import FilterTimeInCycle
            check = FilterTimeInCycle(tref=tref)
            found = False
            # search from last cell in the past until root
            cidup = self.colony.rsearch(self.idseq[-1])
            # note that self.idseq[-1] is a colony leaf
            for cid in cidup:
                cell = self.colony.get_node(cid)
                if check(cell):
                    found = True
                    genref = self.colony.level(cid)
                    break
            if not found:
                msg = 'No cell found in ancestry at {}'.format(tref)
                raise NoAncestry(msg)
        return gens - genref

    def get_boolean_tests(self, cset=[]):
        """Set the dictionary of conditions for a given cset

        Parameters
        ----------
        cset : sequence of :class:`FilterSet` instances (default [])

        Returns
        -------
        select_ids : dict
            keys: 'master' + each of cset item string representation (repr),
            values: sequence of booleans, where index in sequence corresponds
                    to index of each cell in self.idseq

        Notes
        -----
        * 'master' is added, where every test is True
        * call this function AFTER timeseries have been built, so that
          filters acting on Observables can work (computation of
          corresponding cell._sdata is being performed in .get_timeseries(obs))
        """
        # initialize select_ids
        # add the master entry (which will be a sequence of True values)
        select_ids = {}
        # master mask gets all True
        select_ids['master'] = np.array([True for _ in self.idseq])
        # add as many entries as there are conditions
        for fset in cset:
            # we have to make a logical AND between different filter types
            col = self.colony
            cont = col.container
            # check True for upstream structures: container, colony, lineage
            boo = (fset.container_filter(cont) and
                   fset.colony_filter(col) and
                   fset.lineage_filter(self))
            # cell selections
            # initialize all to False
            arrbool = np.array(len(self.idseq) * [False])
            # perform tests only if upstream tests where True
            if boo:
                for index, cid in enumerate(self.idseq):
                    cell = col.get_node(cid)
                    arrbool[index] = fset.cell_filter(cell)
            select_ids[repr(fset)] = arrbool
        return select_ids

    def get_timeseries(self, obs, cset=[]):
        """Contructs timeseries.

        Parameters
        ----------
        obs : :class:`Observable` instance
        cset: sequence of :class:`FilterSet` instances (default [])

        Returns
        -------
        TimeSeries instance
        """
        # check for supplementary observables to be computed
        suppl_obs = []
        for filt in cset:
            # simple filter: ceck for hidden _obs attributes
            if hasattr(filt, '_obs'):
                if isinstance(filt._obs, Observable):
                    suppl_obs.append(filt._obs)
                elif isinstance(filt._obs, collections.Iterable):
                    for item in filt._obs:
                        if isinstance(item, Observable):
                            suppl_obs.append(item)
        # compute suppl obs for all cells in lineage
        if suppl_obs:
            for cid in self.idseq:
                cell = self.colony.get_node(cid)
                for sobs in suppl_obs:
                    _ = cell.build(sobs)

        # build timeseries depending on obs mode
        if obs.mode == 'dynamics':
            return self.get_continuous_timeseries(obs, cset)
        else:
            return self.get_cyclized_timeseries(obs, cset)

    def get_cyclized_timeseries(self, obs, cset=[]):
        """Constructs timeseries for cell cycle observables.

        Parameters
        ----------
        obs : :class:`Observable` instance
            mode should be different from 'dynamics'
        cset: sequence of :class:`FilterSet` instances (default [])

        Returns
        -------
        :class:`TimeSeries` instance
        """
        label = obs.label()
        if obs.timing == 'g':
            try:
                gens = self.get_generations(tref=obs.tref)
            except NoAncestry:
                # return empty TimeSeries
                # TODO : should be filtered upstream?
                new = TimeSeries(label=label, ts=[], ids=self.idseq[:],
                                 index_cycles=[None for _ in self.idseq],
                                 select_ids=self.get_boolean_tests(cset))
                return new

        # browse ids and retrieve stuff
        index_cycles = []
        cts = []
        time_bounds = []
        for index, cid in enumerate(self.idseq):
            cell = self.colony.get_node(cid)
            if cell.birth_time is not None:
                tleft = cell.birth_time
            elif len(cell.data) > 0:
                tleft = np.amin(cell.data['time'])
            else:
                tleft = np.infty
            if cell.division_time is not None:
                tright = cell.division_time
            elif len(cell.data) > 0:
                tright = np.amax(cell.data['time'])
            else:
                tright = - np.infty
            time_bounds.append((tleft, tright))
            value = cell.build(obs)
            # time value
            if obs.timing == 'b':
                tt = cell.birth_time
            elif obs.timing == 'd':
                tt = cell.division_time
            elif obs.timing == 'm':
                try:
                    tt = (cell.division_time + cell.birth_time)/2.
                except TypeError:
                    tt = None
            elif obs.timing == 'g':
                tt = gens[index]
            # append new value
            cts.append((tt, value, int(cell.identifier)))
            index_cycles.append((index, index))

        select_ids = self.get_boolean_tests(cset)
        cts = np.array(cts, dtype=[('time', 'f8'), (label, 'f8'),
                                   ('cellID', 'u2')])
        new = TimeSeries(label=label, ts=cts, ids=self.idseq[:],
                         time_bounds=time_bounds,
                         index_cycles=index_cycles, select_ids=select_ids)
        return new

    def get_continuous_timeseries(self, obs, cset=[]):
        """Method that computes timeseries associated to observable obs.

        Parameters
        ----------
        obs : :class:`Observable` instance
            It defines how to get observable from Numpy structured arrays
        cset: sequence of :class:`FilterSet` instances (default [])

        Returns
        -------
        :class:`TimeSeries` instance
        """
        # build timeseries by concatenating each cell's timeseries
        # when time_window is called, some values from the estimate at cell <c>
        # are in fact evaluated in its parent cell time range
        ts = []
        arrays = []
        valid_previous_cell = False
        # lineage conditions mask
        time_bounds = []
        for cid in self.idseq:
            cell = self.colony.get_node(cid)
            if cell.birth_time is not None:
                tleft = cell.birth_time
            elif len(cell.data) > 0:
                tleft = np.amin(cell.data['time'])
            else:
                tleft = np.infty
            if cell.division_time is not None:
                tright = cell.division_time
            elif len(cell.data) > 0:
                tright = np.amax(cell.data['time'])
            else:
                tright = - np.infty
            time_bounds.append((tleft, tright))
            local = cell.build(obs)  # get local timeseries
            # cut extrapolated values if there is no previous cell
            if len(local) > 0:
                # this is to dismiss data used in other lineages
                if (not valid_previous_cell) and (cell.birth_time is not None):
                    for index, time in enumerate(local['time']):
                        if time >= cell.birth_time:
                            break
                    local = local[index:]
            arrays.append(local)
            # id is added even if not data is reported BY THIS CELL
            # but its daughter cell potentially can report for data
            # in THIS CELL time range due to the local_build() procedure
            valid_previous_cell = cell.data is not None and len(cell.data) > 0

        # try to identify closest frames to cell birth/division
        index_cycles = []

        ts = np.concatenate(arrays)
        if len(ts) > 0:

            # get array of times
            times = ts['time']
            frames = len(times)
            index = 0

            # report indices of timeseries for cell's time range
            index_ts_birth = None
            index_ts_division = None
            for cid in self.idseq:
                cell = self.colony.get_node(cid)
                # cell without data: no range for
                if cell.data is None or len(cell.data) == 0:
                    index_cycles.append(None)
                    continue  # move to next cell

                # last cell may not have data in timeseries
                if cell.birth_time is not None:
                    if cell.birth_time > times[-1]:
                        index_cycles.append(None)
                        continue

                    # if cell.birth_time is defined, look for first frame
                    while index < frames and times[index] < cell.birth_time:
                        index += 1
                    if index < frames:
                        index_ts_birth = index  # index of cell's first frame

                # if no birth is reported, then index_ts_birth is kept to None
                # just checking that cell's data is not beyond time values
                elif cell.data['time'][0] > times[-1]:
                    index_cycles.append(None)
                    continue

                # if cell.division_time is defined, look for last frame
                if cell.division_time is not None:
                    while index < frames and times[index] < cell.division_time:
                        index += 1
                    if index <= frames:
                        index_ts_division = index - 1  # index of cell's last frame
                # if it's not, it means it is last cell in lineage
                else:
                    index_ts_division = None  # slice til the end
                index_cycles.append((index_ts_birth, index_ts_division))
        # if no data has been retrieved, there might still be cells in the lineage
        # but with too few points to
        else:
            index_cycles = [None for cid in self.idseq]

        select_ids = self.get_boolean_tests(cset)

        # return a TimeSeries instance
        timeseries = TimeSeries(label=obs.label(),
                                ts=ts,
                                ids=self.idseq[:],
                                time_bounds=time_bounds,
                                index_cycles=index_cycles,
                                select_ids=select_ids)

        return timeseries

    def split(self):
        "Returns slices to retrieve each cell's data"
        slices = []
        start = 0
        refid = self.data[0]['cellID']
        for index, line in enumerate(self.data[1:], start=1):
            cid = line['cellID']
            if cid != refid:
                stop = index
                slices.append(slice(start, stop))
                start = index
                refid = cid
        stop = None
        slices.append(slice(start, stop))
        return slices

    def __repr__(self):
        label = ''
        if self.colony is not None:
            if self.colony.container is not None:
                label += 'From Container: {}'.format(self.colony.container.label)
            label += ', colony with root id: {}'.format(self.colony.root)
        label += '\nCell ids: {}'.format(' > '.join(self.idseq)) + '\n'
        return label


def get_division_timing(idseq, tree):
    timings = []
    for cid in idseq:
        node = tree.get_node(cid)
        timings.append(node.division_time)
    return timings


def get_lineage_data(idseq, tree, subcellcycle_smooth=False,
                     subcellcycle_smooth_timescale=20.):
    """Builds a numpy array that concatenate data along lineage of a given tree

    Arguments
    ---------
    idline -- list of node identifiers
    tree -- treelib.Tree instance

    Returns
    -------
    ndarray (structured array)
    """
    data = []
    for cid in idseq:
        d = tree.get_node(cid).data
        if d is not None:
            data.append(d)
    if data:
        return np.concatenate(data)
    else:
        return None
