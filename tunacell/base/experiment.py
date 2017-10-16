#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module implements the first core class: Experiment,
and functions to parse containers, retrieve and build data.

Each experiment consists of multiple containers where
data has been stored under container folders. A container may
correspond to a single field of view, to a subset thereof (e.g. a single
channel).

Such container must meet two requirements:
    1. Cell identifiers are unique within this container;
    2. Lineage reconstruction is defined within a single container.

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

from tunacell.base.container import Container
from tunacell.filters.main import FilterSet
from tunacell.io import text, metadata


class ParsingExperimentError(Exception):
    pass


class FiletypeError(ParsingExperimentError):
    pass


class MissingFile():
    pass


class Experiment(object):
    """General class that stores experiment details.

    Creates an Experiment instance from reading a file, records path, filetype,
    reads metadata, stores the list of container containers.

    Parameters
    ----------
    path : str
        path to experiment root file
    filetype -- str {None, 'text', 'h5'}
        leave to None for automatic detection.

    Attributes
    ----------
    abspath : str
        absolute path on disk of main directory for text containers
    label : str
        experiment label
    filetype : str {'text', 'simu'}
        one of the available file type ('simu' is not a filetype per se...)
    fset : :class:`FilterSet` instance
        filterset to be applied when parsing data
    datatype : Numpy.dtype instance
        provides the datatype of raw data stored in each Cell instance .data
        This attribute is defined only for text filetype, when a descriptor
        file is associated to the experiment.
    metadata: Metadata instance
        experiment metadata
    containers : list of str
        sequence of container filenames for text containers
    period: float
        time interval between two successive aquisitions (this should be
        defined in the experiment metadata)

    Methods
    -------
    iter_containers(self, read=True, build=True, prefilt=None,
                    extend_observables=False, report_NaNs=True,
                    size=None, shuffle=False)
        main API to browse containers. MUST BE DEFINED FOR USE IN HIGH LEVEL
        PARSER API CLASS.
    """

    def __init__(self, path='.', filetype=None, filter_set=None):
        self.abspath = None
        self.label = None
        self.datatype = None  # Will be updated for text filetype
        self._containers = []
        self.metadata = None
        self.period = None
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
            self.load_from_text()
        else:
            raise FiletypeError('Filetype not recognized')
        if isinstance(filter_set, FilterSet):
            self._fset = filter_set
        else:
            self._fset = FilterSet()  # default filter: all TRUE
        return

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
        return
    
    @property
    def analysis_path(self):
        """Get analysis path (with appropriate filterset path)"""
        analysis_path = text.get_analysis_path(self, write=True)
        index, filterset_path = text.get_filter_path(analysis_path, self.fset, write=True)
        return filterset_path

    def load_from_text(self):
        """Load parameters from text filetype"""
        # get list of container files
        fns = text.container_filename_parser(self.abspath)
        bns = [os.path.splitext(bn)[0] for bn in fns]
        # extract only container labels
        self.containers = bns
        # get metadata
        fn = text._check_up('metadata.csv', self.abspath, level=2)
        df = metadata.load_from_csv(fn, sep=',')
        meta = metadata.fill_rows(df, self.label, self.containers)
        self.metadata = meta
        self.period = metadata.get_period(meta, self.label)
        descriptor_file = text._check_up('descriptor.csv', self.abspath, 2)
        datatype = text.datatype_parser(descriptor_file)
        self.datatype = datatype
        return

    @property
    def containers(self):
        return self._containers

    @containers.setter
    def containers(self, value):
        if isinstance(value, list):
            self._containers = value
        elif isinstance(value, str):
            self._containers.append(value)
        return

    def info(self):
        """Show informations about experiment"""
        msg = 'Experiment root: {}\n'.format(self.abspath)
        msg += '({} containers)\n'.format(len(self.containers))
        msg += 'Filetype: {}\n'.format(self.filetype)
        msg += repr(self.metadata.loc[self.label])
        return msg

    def __repr__(self):
        msg = 'Experiment root: {}'.format(self.abspath)
        msg += '\nContainers:\n'
        count = 0
        for fn in self.containers:
            count += 1
            if count > 10:
                msg += '\t...\n'
                break
            msg += '\t' + fn + '\n'
        msg += '({} containers)\n'.format(len(self.containers))
        msg += 'Filetype: {}\n'.format(self.filetype)
        msg += repr(self.metadata.loc[self.label])
        return msg

    def iter_containers(self, read=True, build=True, prefilt=None,
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
        prefilt : FilterCell instance (default None)
        apply_container_filter : bool (default True)
            whether to apply container filter defined in filter_set self.fset
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
        if prefilt is None:
            init_cell_filter = self.fset.cell_filter
        else:
            init_cell_filter = prefilt
        containers = self.containers[:]
        if shuffle:
            random.shuffle(containers)
        if size is None:
            size = len(self.containers)
        for index, label in enumerate(containers, start=1):
            if index > size:
                break
            container = Container(label, exp=self)
            if read:
                container.read_data(build=build, prefilt=init_cell_filter,
                                    extend_observables=extend_observables,
                                    report_NaNs=report_NaNs)
            # apply filtering on containers
            if apply_container_filter and self.fset.container_filter(container):
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
        if self.filetype == 'text':
            found = False
            for item in self.containers:
                if label == item:
                    found = True
                    break
            if found:
                container = Container(label, exp=self)
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

    def iter_colonies(self, size=None, shuffle=False):
        """Iterate through valid colonies.

        Parameters
        ----------
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
        if size is not None:
            count = 0  # count colonies
            for container in self.iter_containers(shuffle=shuffle):
                for colony in container.iter_colonies(filt=colfilt,
                                                      shuffle=shuffle):
                    yield colony
                    count += 1
                    if count >= size:
                        break
                if count >= size:
                    break
        else:
            for container in self.iter_containers(shuffle=shuffle):
                for colony in container.iter_colonies(filt=colfilt,
                                                      shuffle=shuffle):
                    yield colony
        return

    def iter_lineages(self, size=None, shuffle=False):
        """Iterate through valid lineages.

        Parameters
        ----------
        size : int (default None)
            limit the number of lineages to size. Works only in mode='all'
        shuffle : bool (default False)
            whether to shuffle the ordering of lineages when mode='all'

        Yields
        ------
        lineage : :class:`Lineage` instance
            filtering removed outlier cells, containers, colonies, and lineages
        """
        lin_filt = self.fset.lineage_filter
        if size is not None:
            count = 0
            for colony in self.iter_colonies(shuffle=shuffle):
                for lineage in colony.iter_lineages(filt=lin_filt,
                                                    shuffle=shuffle):
                    yield lineage
                    count += 1
                    if count >= size:
                        break
                if count >= size:
                    break
        else:
            for colony in self.iter_colonies(shuffle=shuffle):
                for lineage in colony.iter_lineages(filt=lin_filt,
                                                    shuffle=shuffle):
                    yield lineage
        return

    def iter_cells(self, size=None, shuffle=False):
        """Iterate through valid cells.

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

    def raw_text_export(self, path='.'):
        """Export raw data as text containers in correct directory structure.

        Parameters
        ----------
        path : str
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
        fn = os.path.join(exp_path, 'metadata.csv')
        self.metadata.to_csv(fn, index=False)
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
