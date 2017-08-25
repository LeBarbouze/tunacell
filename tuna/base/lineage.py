#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module defines how to handle lineages.
"""
from __future__ import print_function

import numpy as np
import random
import collections
from tuna.datatools import Coordinates

from tuna.observable import Observable, FunctionalObservable

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
    get_boolean_tests(cset=[])
        Performs boolean test over cells in the sequence
    get_timeseries(obs, cset=[])
        Returns :class:`Timeseries` instance of :class:`Observable`
        :param:`obs` in current lineage

    """

    def __init__(self, tree, identifier_sequence):
        self.colony = tree
        self.idseq = identifier_sequence
        self.cellseq = [tree.get_node(cid) for cid in self.idseq]  # loops
        self.division_timings = get_division_timing(self.idseq, tree)
        return

    # DEPRECATED: CHECK WHERE IT IS USED
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
            arrbool = np.array(len(self.idseq) * [False, ])
            # perform tests only if upstream tests where True
            if boo:
                for index, cell in enumerate(self.cellseq):
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
        label = obs.label
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
        # check for supplentary observables when obs is FunctionalObservable
        if isinstance(obs, FunctionalObservable):
            suppl_obs.extend(obs.observables)
        elif isinstance(obs, Observable):
            suppl_obs.append(obs)

        # compute timelapsed obs for all cells in lineage
        for cell in self.cellseq:
            for sobs in suppl_obs:
                cell.build(sobs.as_timelapse())
        # now that all timelapse observables have been computed,
        # compute those that are of cell-cycle mode
        for cell in self.cellseq:
            for sobs in suppl_obs:
                if sobs.mode != 'dynamics':
                    cell.compute_cyclized(sobs)

        time_bounds = []
        for cell in self.cellseq:
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

        # perform functional operation
        if isinstance(obs, FunctionalObservable):
            # define new sdata with functional form
            for cell in self.cellseq:
                arrays = [cell._sdata[item.label] for item in obs.observables]
                result_array = obs.f(*arrays)
                cell._sdata[obs.label] = result_array

        # boolean tests
        select_ids = self.get_boolean_tests(cset)

        arrays = []
        index_cycles = []
        # at this point all _sdata are ready for action. Distinguish modes
        if obs.mode == 'dynamics':
            # build array
            count = 0
            for cell in self.cellseq:
                if len(cell.data) > 0:
                    coords = Coordinates(cell.data['time'], cell._sdata[label])
                    arrays.append(coords.as_array(x_name='time', y_name=label))
                    size = len(arrays[-1])
                    index_cycles.append((count, count + size - 1))
                    count += size
                else:
                    index_cycles.append(None)
            ts = np.concatenate(arrays)
            # return a TimeSeries instance
            timeseries = TimeSeries(label=label,
                                    ts=ts,
                                    ids=self.idseq[:],
                                    time_bounds=time_bounds,
                                    index_cycles=index_cycles,
                                    select_ids=select_ids)
        # otherwise it's of 'cycle' mode
        else:
            for index, cell in enumerate(self.cellseq):
                if timing == 'g':
                    try:
                        gens = self.get_generations(tref=obs.tref)
                        tt = gens[index]
                    except NoAncestry:
                        # return empty TimeSeries
                        # TODO : should be filtered upstream?
                        new = TimeSeries(label=label, ts=[], ids=self.idseq[:],
                                         index_cycles=[None for _ in self.idseq],
                                         select_ids=select_ids)
                        return new
                # time value
                elif timing == 'b':
                    tt = cell.birth_time
                elif timing == 'd':
                    tt = cell.division_time
                elif timing == 'm':
                    try:
                        tt = (cell.division_time + cell.birth_time)/2.
                    except TypeError:
                        tt = None

                # append new value
                arrays.append((tt, cell._sdata[label]))
                index_cycles.append((index, index))

            ts = np.array(arrays, dtype=[('time', 'f8'), (label, 'f8')])
            timeseries = TimeSeries(label=label, ts=ts, ids=self.idseq[:],
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
