#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module provides elements for parsing data manually, i.e. getting
a handful list of samples and extract specific structures (colony, lineage,
cell) from such samples.

* :class:`Parser`: handles how to parse an experiment with a given filterset

"""
from __future__ import print_function

import os
import numpy as np
import warnings

from tunacell.filters.main import FilterSet

from tunacell.base.experiment import Experiment
from tunacell.base.container import ParsingContainerError
from tunacell.base.lineage import Lineage


class Parser(object):
    """Defines how user wants to parse data.

    Parameters
    ----------
    exp : :class:`Experiment` instance
    filter_set : :class:`FilterSet` instance
        this is the set of filters used to read/build data, used for
        for iterators
        (usually, only .cell_filter and .container_filter are used)
    """

    def __init__(self, exp=None, filter_set=None):
        if exp is None:
            print('Use parser.load_experiment() to load from path to file')
        else:
            if isinstance(exp, Experiment):
                self._experiment = exp
            else:
                self.load_experiment(exp)
        if isinstance(filter_set, FilterSet):
            self.fset = filter_set
        elif filter_set is None:
            fset = FilterSet()
            self.fset = fset
        else:
            raise TypeError('filter_set is not a FilterSet instance')
        self._sample_list = []  # used for small sample visualization
        return

    # GETTING/SETTING EXPERIMENT
    @property
    def experiment(self):
        return self._experiment

    @experiment.setter
    def experiment(self, exp):
        """Set experiment.

        Parameters
        ----------
        exp : Experiment instance
        """
        self._experiment = exp
        return

    def load_experiment(self, path, filetype='text'):
        """Loads an experiment from path to file.

        Parameters
        ----------
        path : str
            path to root directory ('text'), or to datafile ('h5')
        filetype : str {'text', 'h5}
        """
        exp = Experiment(path=path, filetype=filetype)
        self.experiment = exp
        return

    def set_filter(self, fset):
        """Set build filter.

        Parameters
        ----------
        fset : :class:`FilterSet` instance
        """
        self.experiment.fset = fset
        return

    # ADDING SAMPLES
    def _add_atomic_sample(self, sample_id):
        """Add atomic sample.

        Parameters
        ----------
        sample_id : dict
            keys: 'container_label', 'cellID'
        """
        if sample_id not in self._sample_list:
            self._sample_list.append(sample_id)
        else:
            print('Sample {} already stored'.format(sample_id))
        return

    def _add_random_sample(self, container_label=None):
        """Add random sample.

        When parameter is set to None, a container is chosen randomly, and
        then a cell is chosen randomly within this container. Otherwise the
        requested container is parsed to draw a random cell out of it.

        Parameters
        ----------
        container_label : str (default None)
            must be a valid container label.
        """
        exp = self.experiment
        filtcell = self.fset.cell_filter
        # 1. container is None
        item = {}
        if container_label is None:
            empty = True
            while empty:
                label = np.random.choice(exp.containers)
                container = exp.get_container(label, read=True,
                                              build=True, prefilt=filtcell,
                                              extend_observables=False,
                                              report_NaNs=True)
                empty = len(container.cells) == 0
        else:
            try:
                err = None
                container = exp.get_container(container_label, read=True,
                                              build=True, prefilt=filtcell,
                                              extend_observables=False,
                                              report_NaNs=True)
                label = container.label
            except ParsingContainerError as e:
                err = e
            if err is not None or len(container.cells) == 0:
                if err is not None:
                    msg = err
                else:
                    msg = "Er, there's no cell in this container, brah*"
                    msg += "\nI'll choose it randomly for you then. Peace*."
                    msg += "\n*If you find this message, all apologies..."
                warnings.warn(msg)
                self._add_random_sample(None)
                return
        cell = np.random.choice(container.cells)
        item['container_label'] = label
        item['cellID'] = cell.identifier
        # we're just checking that current item has not be drawn previously
        if item not in self._sample_list:
            self._add_atomic_sample(item)
        # try again dude
        # COMMENTED SINCE IT CAN LOOP INDEFINITELY
#        else:
#            self._add_random_sample(container_label=container_label)
        return

    def add_sample(self, *args):
        """Add sample to sample list.

        Parameters
        ----------
        args : list
            list of items such as: integer, strings, couple, and/or dict
            An integer denotes the number of sample_ids to be chosen randomly
            (if many integers are given, only the first is used)
            A string will result in loading the corresponding Container,
            with a cell identifier randomly chosen.
            A couple (str, cellID) denotes (container_label, cell_identifier)
            A dictionary should provide 'container' key, and 'cellID' key
        """
        exp = self.experiment
        filtcell = self.fset.cell_filter
        item = {}
        for arg in args:
            if isinstance(arg, int):
                for count in range(arg):
                    self._add_random_sample()
            # argument is container label
            elif isinstance(arg, str):
                self._add_random_sample(container_label=arg)
            # argument is a tuple or a dict already
            else:
                if isinstance(arg, tuple):
                    # try opening container and find cellID
                    label, cid = arg[:2]
                    cid = str(cid)  # just in case
                    label, ext = os.path.splitext(os.path.basename(label))
                elif isinstance(arg, dict):
                    try:
                        label = arg['container_label']
                        cid = str(arg['cellID'])
                    except KeyError as k:
                        msg = k
                        msg += '\nCheck input.'
                        warnings.warn(msg)
                        return
                label = label.replace('data_', '')
                try:
                    # a Parsing error is thrown if label is incorrect
                    container = exp.get_container(label, read=True,
                                                  build=True, prefilt=filtcell,
                                                  extend_observables=False,
                                                  report_NaNs=True)
                    item['container_label'] = label
                    # we'll throw another Parsing error if cell not found
                    found = False
                    for cell in container.cells:
                        if cid == cell.identifier:
                            found = True
                            break
                    if not found:
                        msg = 'cell {} not found in'.format(cid)
                        msg += ' container {}'.format(label)
                        raise ParsingContainerError(msg)
                    item['cellID'] = cid
                # catching Parsing errors and stopping
                except ParsingContainerError as pe:
                    msg = pe
                    msg += '\nCheck input.'
                    warnings.warn(msg)
                    return  # exit
                self._add_atomic_sample(item)
        return

    @property
    def samples(self):
        return self._sample_list

    def get_sample(self, index, level='cell'):
        """Return sample corresponding to index.

        Parameters
        ----------
        index : int
            index of sample id in self.samples
        level : str {'cell'|'lineage'|'colony'}

        Returns
        -------
        out
            structure level corresponding to sample id
        """
        sample = self.samples[index]
        if level == 'cell':
            out = self.get_cell(sample)
        elif level == 'lineage':
            out = self.get_lineage(sample)
        elif level == 'colony':
            out = self.get_colony(sample)
        return out

    def get_cell(self, sample_id):
        """Get :class:`Cell` instance corresponding to sample_id.

        Parameters
        ----------
        sample_id : dict
            element of self.samples

        Returns
        -------
        cell : :class:`Cell` instance
            corresponding to sample_id
        """
        if isinstance(sample_id, int):
            sample_id = self.samples[sample_id]
        exp = self.experiment
        filtcells = self.fset.cell_filter
        label = sample_id['container_label']
        container = exp.get_container(label, read=True, build=True,
                                      prefilt=filtcells,
                                      extend_observables=True,
                                      report_NaNs=True)
        out = None
        for cell in container.cells:
            if cell.identifier == sample_id['cellID']:
                out = cell
                break
        return out

    def get_colony(self, sample_id):
        """Get :class:`Colony` instance corresponding to sample_id.

        Parameters
        ----------
        sample_id : dict
            element of self.samples

        Returns
        -------
        colony : :class:`Colony` instance
            corresponding to sample_id
        """
        if isinstance(sample_id, int):
            sample_id = self.samples[sample_id]
        exp = self.experiment
        filtcells = self.fset.cell_filter
        label = sample_id['container_label']
        container = exp.get_container(label, read=True, build=True,
                                      prefilt=filtcells,
                                      extend_observables=True,
                                      report_NaNs=True)
        out = None
        for colony in container.trees:
            if colony.contains(sample_id['cellID']):
                out = colony
                break
        return out

    def get_lineage(self, sample_id):
        """Get :class:`Lineage` instance corresponding to sample_id.

        Parameters
        ----------
        sample_id : dict
            element of self.samples

        Returns
        -------
        lineage : :class:` Lineage` instance
            corresponding to sample_id
        """
        if isinstance(sample_id, int):
            sample_id = self.samples[sample_id]
        colony = self.get_colony(sample_id)
        ptl = colony.paths_to_leaves()
        candidates = []
        for path in ptl:
            if sample_id['cellID'] in path:
                candidates.append(path)
        candidates.sort(key=lambda x: len(x), reverse=True)
        lineage = Lineage(colony, candidates[0])
        return lineage

    def get_container(self, sample_id):
        """Get :class:`Container` instance corresponding to sample_id.

        Parameters
        ----------
        sample_id: dict
            element of self.samples

        Returns
        -------
        container : :class:`Container` instance
        """
        if isinstance(sample_id, int):
            sample_id = self.samples[sample_id]
        exp = self.experiment
        filtcells = self.fset.cell_filter
        label = sample_id['container_label']
        container = exp.get_container(label, read=True, build=True,
                                      prefilt=filtcells,
                                      extend_observables=True,
                                      report_NaNs=True)
        return container

    def clear_samples(self):
        """Erase all samples."""
        self._sample_list = []
        return

    def remove_sample(self, index, verbose=True):
        """Remove sample of index in sample list."""
        item = self._sample_list.pop(index)
        if verbose:
            print('Item pointed by: ')
            print('Container: {}, cellID: {} '.format(item['container_label'],
                                                      item['cellID']))
            print('has been removed from sample list.')
        return

    def info_samples(self):
        """Table output showing stored samples."""
        column_size = 12
        msg = ''
        if not self._sample_list:
            msg += 'No samples have been added yet. Use .add_sample().'
        else:
            def formatting(items):
                out = ' | '.join(['{}'.format(item).center(column_size)
                                  for item in items])
                return out
            msg += '\n' + formatting(['index', 'container', 'cell'])
            vsep = ''.join(['-' for space in range(column_size)])
            msg += '\n' + formatting([vsep, vsep, vsep])
            for index, sample_id in enumerate(self._sample_list):
                msg += '\n' + formatting([index,
                                          sample_id['container_label'],
                                          sample_id['cellID']])
        return msg

    def __repr__(self):
        return self.info_samples()


    def iter_containers(self, mode='all', size=None, shuffle=False):
        """Iterate through valid containers.

        .. note:: Deprecated in tunacell 0.0.7
                  mode 'all' is similar to run iter_containers from Experiment
                  instance. This mode will be removed in the future.

        If mode 'all' is chosen, then the iterator browses all files,
        up to the size limit. If mode 'samples' is chosen, then only sample ids
        already stored in Parser instance are browsed. Only valid Container
        instance are yielded (valid under filtering against containers).

        Parameters
        ----------
        mode : str {'all', 'samples'} (default 'all')
            whether to iterate over all possible containers, or restrict to
            containers pointed by parser.samples
        size : int (default None), number of containers to be parsed
        shuffle : bool (default False)
            whether to randomize ordering of containers

        Yields
        ------
        container : :class:`tunacell.core.Container` instance
            filtering removed outlier cells, containers
        """
        exp = self.experiment
        if mode == 'all':
            for container in exp.iter_containers(read=True, build=True,
                                                 prefilt=self.fset.cell_filter,
                                                 extend_observables=True,
                                                 report_NaNs=True,
                                                 size=size,
                                                 shuffle=shuffle):
                if self.fset.container_filter(container):
                    yield container
        elif mode == 'samples':
            count = 0
            for sample_id in self.samples:
                container = self.get_container(sample_id)
                if self.fset.container_filter(container):
                    yield container
                    count += 1
                    if size is not None and count >= size:
                        break
        return

    def iter_colonies(self, mode='all', size=None, shuffle=False):
        """Iterate through valid colonies.

        .. note:: Deprecated in tunacell 0.0.7
                  mode 'all' is similar to run iter_containers from Experiment
                  instance. This mode will be removed in the future.

        Parameters
        ----------
        mode : str {'all', 'samples'} (default 'all')
            whether to iterate over all colonies (up to number limitation), or
            over registered samples
        size : int (default None)
            limit the number of colonies to size. Works only in mode='all'
        shuffle : bool (default False)
            whether to shuffle the ordering of colonies when mode='all'

        Yields
        ------
        colony : :class:`Colony` instance
            filtering removed outlier cells, containers, and colonies
        """
        colfilt = self.fset.colony_filter
        if mode == 'all':
            if size is not None:
                count = 0  # count colonies
                for container in self.iter_containers(mode='all',
                                                      shuffle=shuffle):
                    for colony in container.iter_colonies(filt=colfilt,
                                                          shuffle=shuffle):
                        yield colony
                        count += 1
                        if count >= size:
                            break
                    if count >= size:
                        break
            else:
                for container in self.iter_containers(mode='all',
                                                      shuffle=shuffle):
                    for colony in container.iter_colonies(filt=colfilt,
                                                          shuffle=shuffle):
                        yield colony
        elif mode == 'samples':
            count = 0
            for sample_id in self.samples:
                colony = self.get_colony(sample_id)
                if colfilt(colony):
                    yield colony
                    count += 1
                    if size is not None and count >= size:
                        break
        return

    def iter_lineages(self, mode='all', size=None, shuffle=False):
        """Iterate through valid lineages.

        .. note:: Deprecated in tunacell 0.0.7
                  mode 'all' is similar to run iter_containers from Experiment
                  instance. This mode will be removed in the future.

        Parameters
        ----------
        mode : str {'all', 'samples'} (default 'all')
            whether to iterate over all lineages (up to number limitation), or
            over registered samples
        size : int (default None)
            limit the number of lineages to size. Works only in mode='all'
        shuffle : bool (default False)
            whether to shuffle the ordering of lineages when mode='all'

        Yields
        ------
        lineage : :class:`Lineage` instance
            filtering removed outlier cells, containers, colonies, and lineages
        """
        if mode == 'all':
            if size is not None:
                count = 0
                for colony in self.iter_colonies(mode='all', shuffle=shuffle):
                    for lineage in colony.iter_lineages(filt=self.fset.lineage_filter,
                                                        shuffle=shuffle):
                        yield lineage
                        count += 1
                        if count >= size:
                            break
                    if count >= size:
                        break
            else:
                for colony in self.iter_colonies(mode='all', shuffle=shuffle):
                    for lineage in colony.iter_lineages(filt=self.fset.lineage_filter,
                                                        shuffle=shuffle):
                        yield lineage
        elif mode == 'samples':
            count = 0
            for sample_id in self.samples:
                lineage = self.get_lineage(sample_id)
                if self.fset.lineage_filter(lineage):
                    yield lineage
                    count += 1
                    if size is not None and count >= size:
                        break
        return

    def iter_cells(self, mode='all', size=None, shuffle=False):
        """Iterate through valid cells.

        .. note:: Deprecated in tunacell 0.0.7
                  mode 'all' is similar to run iter_containers from Experiment
                  instance. This mode will be removed in the future.

        Parameters
        ----------
        mode : str {'all', 'samples'} (default 'all')
            whether to iterate over all cells (up to number limitation), or
            over registered samples
        size : int (default None)
            limit the number of lineages to size. Works only in mode='all'
        shuffle : bool (default False)
            whether to shuffle the ordering of lineages when mode='all'

        Yields
        ------
        cell : :class:`Cell` instance
            filtering removed outlier cells, containers, colonies, and lineages
        """
        if mode == 'all':
            if size is not None:
                count = 0
                for lin in self.iter_lineages(mode='all', shuffle=shuffle):
                    for cell in lin.iter_cells(shuffle=shuffle):
                        yield cell
                        count += 1
                        if count >= size:
                            break
                    if count >= size:
                        break
            else:
                for lin in self.iter_lineages(mode='all', shuffle=shuffle):
                    for cell in lin.iter_cells(shuffle=shuffle):
                        yield cell
        elif mode == 'samples':
            count = 0
            for sample_id in self.samples:
                cell = self.get_cell(sample_id)
                yield cell
                count += 1
                if size is not None and count >= size:
                    break
        return
