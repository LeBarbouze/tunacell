#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module implements the first core class: Experiment,
and functions to parse containers, retrieve and build data.

Each experiment consists of multiple containers where
data is stored under container folders. A container may
correspond to a single field of view, to a subset thereof (e.g. a single
channel in microfluidic experiments).

Such containers must meet two requirements:
    1. Cell identifiers are unique within a container;
    2. Lineage reconstruction is defined and performed within a single
       container.

This module stores classes and functions that allow to:
    * explore data structure (and metadata if provided)
    * keep track of every container where to look for data
    * extract data in a Container instance from text containers
    * build cells filiation, store time-lapse microscopy data, build trees
"""
from __future__ import print_function

import os
import random
import warnings
import shutil

import sys
if sys.version_info[0] < 3:
    import pathlib2 as pathlib
else:
    import pathlib

from tqdm import tqdm

from tunacell.base.container import Container
from tunacell.filters.main import FilterSet, FilterTRUE
from tunacell.filters.cells import FilterCell
from tunacell.filters.containers import FilterContainer
from tunacell.filters.trees import FilterTree
from tunacell.filters.lineages import FilterLineage
from tunacell.io import text, supersegger, analysis


class ParsingExperimentError(Exception):
    pass


class FiletypeError(ParsingExperimentError):
    pass


class MissingFile():
    pass


class Experiment(object):
    """General class that stores experiment details.

    Creates an Experiment instance from reading a file, records path, filetype,
    reads metadata, stores the list of containers.

    Parameters
    ----------
    path : str
        path to experiment root file
    filetype -- str {None, 'text', 'supersegger'}
        leave to None for automatic detection.

    Attributes
    ----------
    abspath : str
        absolute path on disk of main directory for text containers
    label : str
        experiment label
    filetype : str {'text', 'supersegger'}
        one of the available file type ('simu' is not a filetype per se...)
    fset : :class:`FilterSet` instance
        filterset to be applied when parsing data
    datatype : Numpy.dtype instance
        provides the datatype of raw data stored in each Cell instance .data
        This attribute is defined only for text filetype, when a descriptor
        file is associated to the experiment.
    metadata: Metadata instance
        experiment metadata
    containers : list of pathlib.Path
        list of absolute paths to containers
    period: float
        time interval between two successive aquisitions (this should be
        defined in the experiment metadata)

    Methods
    -------
    iter_containers(self, read=True, build=True, prefilt=None,
                    extend_observables=False, report_NaNs=True,
                    size=None, shuffle=False)
        browse containers

    """

    def __init__(self, path='.', filetype=None, filter_set=None,
                 count_items=False):
        self.abspath = None
        self.label = None
        self.datatype = None  # Will be updated for text filetype
        self._containers = []
        self._counts = None
        self.metadata = None
        # let's go
        self.abspath = os.path.abspath(os.path.expanduser(path))
        # remove extension
        label, extension = os.path.splitext(os.path.basename(self.abspath))
        self.label = label
        # recognize filetype from file extension
        if filetype is None:
            if extension == '':
                self.filetype = 'text'
        # if it's not recognized, go default
        else:
            self.filetype = filetype
        # different initialization depending on filetype
        if self.filetype == 'text':
            containers = text.find_containers(self.abspath)  # full path
#            basenames = [item.stem for item in containers]  # only labels
            self.containers = containers
            self.metadata = text.find_metadata(self.abspath)
            self.datatype = text.find_datatype(self.abspath)
        elif self.filetype == 'supersegger':
            containers = supersegger.find_containers(self.abspath)
#            basenames = [item.stem for item in containers]
            self.containers = containers
            self.metadata = text.find_metadata(self.abspath)
        else:
            raise FiletypeError('Filetype not recognized')
        self.fset = filter_set
        if count_items:
            self.count_items()

    @property
    def period(self):
        """Return the experimental level period

        The experimental level period is defined as the smallest acquisition
        period over all containers.
        """
        return self.metadata.period

    @property
    def fset(self):
        """Get current FilterSet"""
        return self._fset

    @fset.setter
    def fset(self, value):
        """Set current FilterSet"""
        if isinstance(value, FilterSet):
            self._fset = value
        elif value is None:
            self._fset = FilterSet()  # default filterset: all TRUE
        else:
            warnings.warn('{} is not a FilterSet'.format(value))

    @property
    def analysis_path(self):
        """Get analysis path (with appropriate filterset path)"""
        analysis_path = analysis.get_analysis_path(self, write=True)
        index, filterset_path = analysis.get_filter_path(analysis_path, self.fset, write=True)
        return filterset_path

    def count_items(self, independent=True, seed=None, read=True):
        """Parse data to count items: cells, colonies, lineages, containers

        Parameters
        ----------
        independent : bool {True, False}
            lineage decomposition parameter
        seed : int, or None
            lineage decomposition parameter
        read : bool {True, False}
            try to read it in analysis folder
        """
        # 1. try to read
        try:
            if read:
                # get analysis path

                analysis_path = analysis.get_analysis_path(self, write=False)
                i, filter_path = analysis.get_filter_path(analysis_path, self.fset,
                                                      write=False)
                counts = analysis.read_count_file(filter_path)
                self._counts = counts
            else:
                raise text.MissingFileError  # mock it to go to exception
        except (text.MissingFileError, text.MissingFolderError, text.CorruptedFileError):
            self._count_items(independent=independent, seed=seed, write=True)
        print(self._count_summary())

    def _count_summary(self):
        msg = ('\nCount summary:\n'
               ' - cells : {}'.format(self._counts['cells']) + '\n'
               ' - lineages : {}'.format(self._counts['lineages']) + '\n'
               ' - colonies : {}'.format(self._counts['colonies']) + '\n'
               ' - containers : {}'.format(self._counts['containers']))
        return msg

    def _count_items(self, independent=True, seed=None, write=True):
        """Hidden method, see above for parameters"""
        counts = count_items(self, independent_decomposition=independent, seed=seed)
        # save the counts
        analysis.write_count_file(self.analysis_path, counts)
        if write:
            self._counts = counts

    def _erase_count_file(self):
        try:
            analysis_path = analysis.get_analysis_path(self, write=False)
            i, filter_path = analysis.get_filter_path(analysis_path, self.fset,
                                                      write=False)
            filename = os.path.join(filter_path, '.counts.yml')
            if os.path.exists(filename):
                os.remove(filename)
        except (text.MissingFileError, text.MissingFolderError):
            pass

    @property
    def containers(self):
        return self._containers

    @containers.setter
    def containers(self, value):
        if isinstance(value, list):
            self._containers = value
        elif isinstance(value, pathlib.Path):
            self._containers.append(value)
        return

    def info(self):
        """Show informations about experiment"""
        msg = 'Experiment root: {}\n'.format(self.abspath)
        msg += '({} containers)\n'.format(len(self.containers))
        msg += 'Filetype: {}\n'.format(self.filetype)
        msg += repr(self.metadata)
        return msg

    def __repr__(self):
        msg = 'Experiment root: {}'.format(self.abspath)
        msg += '\nContainers:\n'
        count = 0
        for fn in self.containers:
            count += 1
            if count > 5:
                msg += '\t...\n'
                break
            msg += '\t' + str(fn) + '\n'
        msg += '\t({} containers)\n'.format(len(self.containers))
#        msg += 'Filetype: {}\n'.format(self.filetype)
        msg += repr(self.metadata)
        if self._counts is not None:
            msg += self._count_summary()
        return msg

    def iter_containers(self, read=True, build=True,
                        filter_for_cells='from_fset',
                        filter_for_containers='from_fset',
                        apply_container_filter=True,
                        extend_observables=False, report_NaNs=True,
                        size=None, shuffle=False):
        """Iterator over containers.

        Parameters
        ----------
        size : int (default None)
            number of containers to be parsed
        read : bool (default True)
            whether to read data and extract Cell instances
        build : bool (default True), called only if `read` is True
            whether to build colonies
        filter_for_cells : FilterCell instance, or str {'from_fset', 'none'}
            filter applied to cells when data files are parsed
        filter_for_containers : FilterContainer instance or str {'from_fset', 'none'}
            filter applied to containers when data files are parsed
        extend_observables : bool (default False)
            whether to construct secondary observables from raw data
        report_NaNs : bool (default True)
            whether to report for NaNs found in data
        shuffle : bool (default False)
            when `size` is set to a number, whether to randomize ordering of
            upcoming containers

        Returns
        -------
        iterator iver Container instances of current Experiment instance.
        """
        if filter_for_cells is 'from_fset':
            cell_filter = self.fset.cell_filter
        elif filter_for_cells is None or filter_for_cells is 'none':
            cell_filter = FilterTRUE()
        elif isinstance(filter_for_cells, FilterCell):
            cell_filter = filter_for_cells
        else:
            raise ValueError('"filter_for_cells" parameter not recognized')
        if filter_for_containers is 'from_fset':
            container_filter = self.fset.container_filter
        elif filter_for_containers is None or filter_for_containers is 'none':
            container_filter = FilterTRUE()
        elif isinstance(filter_for_containers, FilterContainer):
            container_filter = filter_for_containers
        containers = self.containers[:]
        if shuffle:
            random.shuffle(containers)
        if size is None:
            size = len(self.containers)
        for index, path in enumerate(containers, start=1):
            if index > size:
                break
            container = Container(path, exp=self)
            if read:
                container.read_data(build=build, prefilt=cell_filter,
                                    extend_observables=extend_observables,
                                    report_NaNs=report_NaNs)
            # apply filtering on containers
            if container_filter(container):
                yield container
        return

    def get_container(self, label,
                      read=True, build=True, prefilt=None,
                      extend_observables=False, report_NaNs=True):
        """Open specified container from this experiment.

        Parameters
        ----------
        label : str
            name of the container file to be opened
        read : bool (default True)
            whether to read data and extract Cell instances list
        build : bool (default True)
            when `read` option is active, whether to build Colony instances
        extend_observables : bool (default False)
            whether to compute secondary observables from raw data
        report_NaNs: bool (default True)
            whether to report for NaNs found in data

        Returns
        -------
        container : Container instance

        Raises
        ------
        ParsingExperimentError : when no container corresponds in this exp
        ParsingContainerError: when despite of existing container filename,
            parsing of container failed and nothing is loaded
        """
        if self.filetype == 'text' or self.filetype == 'supersegger':
            found = False
            for path in self.containers:
                if label == path.stem:
                    found = True
                    break
            if found:
                container = Container(path, exp=self)
            else:
                msg = 'Filename error: {}'.format(label)
                msg += ' does not correspond to any container file.'
                raise ParsingExperimentError(msg)
        if container is not None:
            if read:
                container.read_data(build=build, prefilt=prefilt,
                                    extend_observables=extend_observables,
                                    report_NaNs=report_NaNs)
            return container
        else:
            msg = 'Filename corresponds to a container'
            msg += ' but somehow container initialization failed'
            raise ParsingExperimentError(msg)

    def iter_colonies(self, filter_for_colonies='from_fset',
                      size=None, shuffle=False):
        """Iterate through valid colonies.

        Parameters
        ----------
        filter_for_colonies : FilterTree instance or str {'from_fset', 'none'}
        size : int (default None)
            limit the number of colonies to size. Works only in mode='all'
        shuffle : bool (default False)
            whether to shuffle the ordering of colonies when mode='all'

        Yields
        ------
        colony : :class:`Colony` instance
            filtering removed outlier cells, containers, and colonies
        """
        if filter_for_colonies is 'from_fset':
            colony_filter = self.fset.colony_filter
        elif filter_for_colonies is None or filter_for_colonies is 'none':
            colony_filter = FilterTRUE()
        elif isinstance(filter_for_colonies, FilterTree):
            colony_filter = filter_for_colonies
        elif isinstance(filter_for_colonies, FilterTRUE):
            colony_filter = filter_for_colonies
        else:
            raise ValueError('"filter_for_colonies" parameter not recognized')
#        colfilt = self.fset.colony_filter
        if size is not None:
            count = 0  # count colonies
            for container in self.iter_containers(shuffle=shuffle):
                for colony in container.iter_colonies(filter_for_colonies=colony_filter,
                                                      shuffle=shuffle):
                    yield colony
                    count += 1
                    if count >= size:
                        break
                if count >= size:
                    break
        else:
            for container in self.iter_containers(shuffle=shuffle):
                for colony in container.iter_colonies(filter_for_colonies=colony_filter,
                                                      shuffle=shuffle):
                    yield colony
        return

    def iter_lineages(self, filter_for_lineages='from_fset',
                      size=None, shuffle=False):
        """Iterate through valid lineages.

        Parameters
        ----------
        filter_for_lineages : FilterLineage instance or str {'from_fset', 'none'}
            filter lineages
        size : int (default None)
            limit the number of lineages to size. Works only in mode='all'
        shuffle : bool (default False)
            whether to shuffle the ordering of lineages when mode='all'

        Yields
        ------
        lineage : :class:`Lineage` instance
            filtering removed outlier cells, containers, colonies, and lineages
        """
        if filter_for_lineages is 'from_fset':
            lineage_filter = self.fset.lineage_filter
        elif filter_for_lineages is None or filter_for_lineages is 'none':
            lineage_filter = FilterTRUE()
        elif isinstance(filter_for_lineages, FilterLineage):
            lineage_filter = filter_for_lineages
        elif isinstance(filter_for_lineages, FilterTRUE):
            lineage_filter = filter_for_lineages
        else:
            raise ValueError('"filter_for_lineages" parameter not recognized')
#        lin_filt = self.fset.lineage_filter
        if size is not None:
            count = 0
            for colony in self.iter_colonies(shuffle=shuffle):
                for lineage in colony.iter_lineages(filter_for_lineages=lineage_filter,
                                                    shuffle=shuffle):
                    yield lineage
                    count += 1
                    if count >= size:
                        break
                if count >= size:
                    break
        else:
            for colony in self.iter_colonies(shuffle=shuffle):
                for lineage in colony.iter_lineages(filter_for_lineages=lineage_filter,
                                                    shuffle=shuffle):
                    yield lineage
        return

    def iter_cells(self, size=None, shuffle=False):
        """Iterate through valid cells.

        Applies all filters defined in fset.

        Parameters
        ----------
        size : int (default None)
            limit the number of lineages to size. Works only in mode='all'
        shuffle : bool (default False)
            whether to shuffle the ordering of lineages when mode='all'

        Yields
        ------
        cell : :class:`Cell` instance
            filtering removed outlier cells, containers, colonies, and lineages
        """
        if size is not None:
            count = 0
            for lin in self.iter_lineages(shuffle=shuffle):
                for cell in lin.iter_cells(shuffle=shuffle):
                    yield cell
                    count += 1
                    if count >= size:
                        break
                if count >= size:
                    break
        else:
            for lin in self.iter_lineages(shuffle=shuffle):
                for cell in lin.iter_cells(shuffle=shuffle):
                    yield cell
        return

    def raw_text_export(self, path='.', metadata_extension='.yml'):
        """Export raw data as text containers in correct directory structure.

        Parameters
        ----------
        path : str
            path to experiment root directory
        metadata_extension : str (default '.yml')
            type of metadata file (now only yaml file works)
        """
        abspath = os.path.abspath(os.path.expanduser(path))
        # check that path is not self.abspath
        if abspath == self.abspath:
            raise ValueError('cannot be the path of original data')
        # create directory structure
        exp_path = os.path.join(abspath, self.label)
        if not os.path.exists(exp_path):
            os.makedirs(exp_path)
        # write metadata file
        if metadata_extension == '.yml':
            fn = os.path.join(exp_path, 'metadata.yml')
            with open(fn, 'w') as f:
                self.metadata.to_yaml(f)
        # write descriptor file
        with open(os.path.join(exp_path, 'descriptor.csv'), 'w') as f:
            for key, value in self.datatype:
                f.write(str(key) + ',' + str(value) + '\n')
        data_path = os.path.join(exp_path, 'containers')
        if not os.path.exists(data_path):
            os.makedirs(data_path)
        else:
            msg = ("Container files already exists.\n"
                   "They will be erased before proceeding.")
            warnings.warn(msg)
            shutil.rmtree(data_path)
            os.makedirs(data_path)
        for container in self.iter_containers():
            container.write_raw_text(data_path)
        return


def count_items(exp, independent_decomposition=True, seed=None):
    """Parse the experiment, with associated FilterSet, and count items"""
    cells = 0
    colonies = 0
    lineages = 0
    containers = 0
    for container in tqdm(exp.iter_containers(read=True, build=True,
                                              filter_for_cells='from_fset',
                                              filter_for_containers='from_fset',
                                              apply_container_filter=True,),
                         total=len(exp.containers), desc='counting items'):
        containers += 1
        for colony in container.iter_colonies(filter_for_colonies='from_fset'):
            colonies += 1
            for lineage in colony.iter_lineages(independent=independent_decomposition, seed=seed):
                lineages += 1
                for cid in lineage.idseq:
                    cells += 1
    counts = {'cells': cells,
              'lineages': lineages,
              'colonies': colonies,
              'containers': containers}
    return counts
