#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module implements the first core class: Experiment,
and functions to parse containers, retrieve and build data.

Data organization
-----------------

Data is stored by experiment. Each experiment consists of multiple containers where
data has been stored. Lowest level containers are containers. A container may
correspond to a single field of view, to a subset thereof (e.g. a single
channel).

Such container must meet two requirements:
    1. Cell identifiers are unique within this container;
    2. Lineage reconstruction is defined within a single container.

This module stores classes and functions that allow to:
    * explore data structure (and metadata if provided)
    * keep track of every container where to look for data
    * extract data in a Container instance
        - from text containers
        - [TODO] from HDF5 storage
    * build cells filiation, store time-lapse microscopy data, build trees
"""
from __future__ import print_function

import os
import numpy as np
import tables
import datetime
import random
import warnings
import shutil

from tuna.base.container import Container

#from tuna.base.metadata import Metadata, get_time_interval
from tuna.io import text, metadata


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
        absolute path on disk of main directory for text containers, or hdf5 file
    label : str
        experiment label
    filetype : str {'text', 'h5', 'simu'}
        one of the available file type ('simu' is not a filetype per se...)
    datatype : Numpy.dtype instance
        provides the datatype of raw data stored in each Cell instance .data
        This attribute is defined only for text filetype, when a descriptor
        file is associated to the experiment.
    metadata: Metadata instance
        experiment metadata
    containers : list of str
        sequence of container filenames for text containers, of container labels for
        h5 containers
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

    def __init__(self,  path='.', filetype=None):
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
            elif extension in ['h5', 'hdf5']:
                self.filetype = 'h5'
        # if it's not recognized, go default
        else:
            self.filetype = filetype
        # different initialization depending on filetype
        if self.filetype == 'text':
            # get list of container files
            fns = text.container_filename_parser(self.abspath)
            # extract only container labels
            self.containers = list(zip(*map(os.path.splitext, fns))[0])

            fn = text._check_up('metadata.csv', self.abspath, level=2)
            df = metadata.load_from_csv(fn, sep=',')
            meta = metadata.fill_rows(df, self.label, self.containers)
            self.metadata = meta
            self.period = metadata.get_period(meta, self.label)
            descriptor_file = text._check_up('descriptor.csv', self.abspath, 2)
            datatype = text.datatype_parser(descriptor_file)
            self.datatype = datatype
        # TODO : review this part when updating HDF5 filetype reading
        elif self.filetype == 'h5':
            pass
            # list tables under the .containers attr for compatibility with text
#            h5file = tables.open_file(self.abspath, mode='r')
#            containers = []
#            for table in h5file.walk_nodes(h5file.root.lineages,
#                                           classname='Table'):
#                containers.append(table._v_name)
#            self.containers = containers
#            # load experiment metadata
#            content = []
#            for label in h5file.root._v_attrs._f_list('user'):
#                content.append((label,
#                                h5file.root._v_attrs.__getitem__(label)))
#            self.metadata = Metadata(content=content)
#            h5file.close()
        else:
            raise FiletypeError('Filetype not recognized')
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
        msg += '\nContainer containers:\n'
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

    def iter_container(self, read=True, build=True, prefilt=None,
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
        containers = self.containers[:]
        if shuffle:
            random.shuffle(containers)
        if size is None:
            size = len(self.containers)
        if self.filetype == 'h5':
            h5file = tables.open_file(self.abspath)
        else:
            h5file = None
        for index, label in enumerate(containers, start=1):
            if index > size:
                break
            container = Container(label, exp=self)
            if read:
                container.read_data(build=build, prefilt=prefilt,
                                    extend_observables=extend_observables,
                                    report_NaNs=report_NaNs)
            yield container
        if h5file is not None:
            h5file.close()
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
        elif self.filetype == 'h5':
            found = False
            for item in self.containers:
                lab = item.replace('data_', '')
                if label == lab:
                    found = True
                    break
            if found:
                container = Container(item, exp=self)
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

    # TODO: update
    def h5_export(self, directory='~', filename=None, overide=False,
                  prefilt=None, testing=True, extend_observables=False,
                  out=False):
        """Export data as a HDF5 archive.

        Parameters
        ----------
        directory -- str, where to store the h5 file
        filename -- str, name of file (without '.h5' extension)
        overide -- bool (default False), overide existing file
        prefilt -- prefiltering function (default None)
        testing -- bool (default True), when True, export only 10 containers
        out -- bool (default False), output or not tables.File object

        Returns
        -------
        if out is True, returns the corresponding tables.File object, that one
        has to close at some point.
        """
        path = os.path.abspath(os.path.expanduser(directory))
        if filename is None:
            bn = os.path.basename(self.abspath)
        else:
            fn, fnext = os.path.splitext(filename)
            bn = os.path.basename(fn)
        fn = os.path.join(path, bn + '.h5')  # h5 filename
        # check whether h5 export already exists
        if os.path.exists(fn) and not overide:
            print('h5 file already exists. Set overide to True to update.')
            return
        # associated log filename
        logfn = os.path.join(path, bn + '.log')
        log = ''
        # associated warnings file
        warnsfn = os.path.join(path, bn + '.warnings')
        warns = ''
        log += 'Log file for h5_export from Experiment instance\n'
        log += repr(self)
        now = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log += '\nExported: {}\n'.format(now)
        log += '\nErrors/warnings are reported in {}\n'.format(warnsfn)
        log += '\nFILTER USED:\n'
        log += repr(prefilt)
        log += '\n'
        # go on
        h5file = tables.open_file(fn, mode='w', title='Experiment file')
        try:
            # add metadata as attributes
            for k, v in self.metadata.content:
                h5file.root._v_attrs.__setattr__(k, v)
            # create the lineages folder where container tables are stored
            lineages = h5file.create_group(h5file.root, 'lineages',
                                           'Microscopy data flat containers')
            # loop over containers
            if testing:
                size = 10  # only 10 containers for testing
                log += '\nTesting mode: only {} containers.\n'.format(size)
            else:
                size = None  # all containers
            log += '\nCOUNTS:\n'
            log += '\n{}\t{}\t{}'.format('Label', 'Cells', 'Trees')
            count_cells, count_trees = 0, 0
            for cont in self.iter_container(size=size, prefilt=prefilt,
                                            extend_observables=extend_observables,
                                            report_NaNs=True):
                if cont.log is not None:
                    warns += cont.log
                arrs = []
                for cell in cont.cells:
                    arrs.append(cell.data)
                arr = np.concatenate(arrs)
                lab = 'Cells from container {}'.format(cont.label)
                tab = h5file.create_table(h5file.root.lineages,
                                          'data_{}'.format(cont.label),
                                          arr.dtype,
                                          lab)
                tab.append(arr)
                tab.flush()
                # logging
                log += '\n{}\t{}\t{}'.format(cont.label, len(cont.cells),
                                             len(cont.trees))
                count_cells += len(cont.cells)
                count_trees += len(cont.trees)
            log += '\n\nTotal\t{}\t{}'.format(count_cells, count_trees)
            # write log file
            with open(logfn, 'w') as f:
                f.write(log)
            # if warnings showed up
            if warns:
                with open(warnsfn, 'w') as f:
                    f.write(warns)
            # close file if no output is needed
            if not out:
                h5file.close()
                return
            else:
                return h5file
        except Exception as e:
            h5file.close()
            raise e

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
        for container in self.iter_container():
            container.write_raw_text(data_path)
        return
