#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module defines how to handle lineages.
"""
from __future__ import print_function

import numpy as np
import random

from tunacell.base.datatools import Coordinates
from tunacell.base.timeseries import TimeSeries
from tunacell.base.observable import Observable, FunctionalObservable


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
            from tunacell.filters.cells import FilterTimeInCycle
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

    def get_timeseries(self, obs, raw_obs=[], func_obs=[], cset=[]):
        """Contructs timeseries.

        Parameters
        ----------
        obs : :class:`Observable` or :class:`FunctionalObservable` instance
            must be an item of raw_obs or an item of func_obs
        raw_obs : list of :class:`Observable` instances
            needed to be computed for filtering or in the case of FunctionalObservable
        func_obs : list of :class:`FunctionalObservable` instances
            needed to be computed for filtering
        cset: sequence of :class:`FilterSet` instances (default [])

        Returns
        -------
        :class:`TimeSeries` instance
            corresponding to obs argument
        """
        label = obs.label  # complicated string
        if obs.name is not None:
            obs_name = obs.name  # simpler string if provided by user
        else:
            obs_name = label
        
        # obs must be either a member of raw_obs, or a member of func_obs
        if isinstance(obs, Observable):
            if obs not in raw_obs:
                raw_obs.append(obs)
        elif isinstance(obs, FunctionalObservable):
            if obs not in func_obs:
                func_obs.append(obs)
        else:
            raise TypeError('obs must be one of {Observable, FunctionalObservable}')

        # compute timelapsed raw obs for all cells in lineage
        for cell in self.cellseq:
            for sobs in raw_obs:
                cell.build(sobs.as_timelapse())
        # now that all timelapse observables have been computed, there cannot
        # be overlap between different cell in data evaluation,
        #and we protect against future build
        time_bounds = []
        for cell in self.cellseq:
            # compute those that are of cell-cycle mode
            for sobs in raw_obs:
                if sobs.mode != 'dynamics':
                    cell.build(sobs)
                # protect against future build for raw observable
                cell.protect_against_build(sobs)
            for fobs in func_obs:
                cell.build(fobs)
                cell.protect_against_build(fobs)

            # collect make time bounds
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

        # boolean tests
        select_ids = self.get_boolean_tests(cset)

        arrays = []
        index_cycles = []
        colony = self.colony
        container = colony.container
        # at this point all _sdata are ready for action. Distinguish modes
        if obs.mode == 'dynamics':
            # get time reference for translating times
            if obs.tref is not None:
                if obs.tref == 'root':
                    root = colony.get_node(colony.root)
                    if root.data is not None and len(root.data) > 0:
                        tref = root.data['time'][-1]  # last time of root
                    else:
                        tref = 0.
                elif isinstance(obs.tref, float) or isinstance(obs.tref, int):
                    tref = float(obs.tref)
                else:
                    tref = 0.
            else:
                tref = 0.
            # build array
            count = 0
            for cell in self.cellseq:
                if len(cell.data) > 0:
                    coords = Coordinates(cell.data['time'] - tref,
                                         cell._sdata[label],
                                         x_name='time',
                                         y_name=obs_name)
                    arrays.append(coords.clear.as_array())  # remove NaNs
                    size = len(arrays[-1])
                    index_cycles.append((count, count + size - 1))
                    count += size
                else:
                    index_cycles.append(None)
            ts = np.concatenate(arrays)
            coords = Coordinates.from_array(ts)
        # otherwise it's of 'cycle' mode
        else:
            for index, cell in enumerate(self.cellseq):
                if obs.timing == 'g':
                    try:
                        gens = self.get_generations(tref=obs.tref)
                        tt = gens[index]
                    except NoAncestry:
                        # return empty TimeSeries
                        # TODO : should be filtered upstream?
                        coords = Coordinates([], [], x_name='g', y_name=obs.name)
                        new = TimeSeries(ts=coords, ids=self.idseq[:],
                                         index_cycles=[None for _ in self.idseq],
                                         select_ids=select_ids,
                                         container_label=container.label,
                                         experiment_label=container.exp.label)
                        return new
                # time value
                elif obs.timing == 'b':
                    tt = cell.birth_time
                elif obs.timing == 'd':
                    tt = cell.division_time
                elif obs.timing == 'm':
                    try:
                        tt = (cell.division_time + cell.birth_time)/2.
                    except TypeError:
                        tt = None

                # append new value
                if tt is None:
                    tt = np.nan
                arrays.append((tt, cell._sdata[label]))
                index_cycles.append((index, index))
            if len(arrays) == 0:
                coords = Coordinates([], [], x_name=obs.timing, y_name=obs.name)
            else:
                coords = Coordinates(*list(zip(*arrays)), x_name=obs.timing, y_name=obs.name)
        timeseries = TimeSeries(ts=coords, ids=self.idseq[:],
                                time_bounds=time_bounds,
                                index_cycles=index_cycles,
                                select_ids=select_ids,
                                container_label=container.label,
                                experiment_label=container.exp.label)
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
