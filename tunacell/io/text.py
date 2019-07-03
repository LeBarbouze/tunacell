#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Modules that defines import from/export to text files
"""
from __future__ import print_function

import os
import glob
import logging
import numpy as np

import sys
if sys.version_info[0] < 3:
    import pathlib2 as pathlib
else:
    import pathlib

from tunacell.io import metadata
from tunacell.base.cell import Cell
from tunacell.base.datatools import compute_secondary_observables



logger = logging.getLogger(__name__)


class TextParsingError(Exception):
    """General class for Errors while parsing text flat files."""
    pass


class MissingFileError(TextParsingError):
    """Class for missing file."""
    pass


class MismatchFileError(TextParsingError):
    """Class for mismatching filename

    Parameters
    ----------
    args : tuple
        arguments passed to Exception class init
    """
    def __init__(self, level='none', *args):
        super().__init__(*args)
        self.level = level


class MissingFolderError(TextParsingError):

    def __init__(self, missing):
        self.missing = missing

    def __str__(self):
        msg = 'Missing folder: {}'.format(self.missing)
        return msg


class CorruptedFileError(TextParsingError):
    """When a file does not contain what it should"""
    pass


class TextFileSystem(TextParsingError):
    """Error when text flat files do not show appropriate structure."""

    def __init__(self, root):
        self.root = root

    def __str__(self):
        label = ''
        label += 'DID NOT FIND APPROPRIATE FILESYSTEM in root\n'
        label += '  > {}\n'.format(self.root)
        label += 'Appropriate structure from root directory:\n'
        label += 'root/\n'
        label += '  descriptor.csv (describe data structure within files)\n'
        label += '  lineages/\n'
        label += '    fov_001.txt\n'
        label += '    fov_002.txt\n'
        label += '    [...]\n'
        return label


def _check_up(filename, path='.', level=2):
    """Check in parent directories (up to level) whether filename exists

    Parameters
    ----------
    filename : str
        name of file
    path : str
        directory to start with (level 0)
    level : int
        how far you'd like to reach upwards (1: parent directory,
        2: parent to parent directory, ...)

    Returns
    -------
    absolute path to file found upward

    Raises
    ------
    :exception:`MissingFileError` when filename has not been found
    """
    rootpath = os.path.abspath(os.path.expanduser(path))
    current = 0
    found = False
    while current <= level:
        current += 1
        path = os.path.join(rootpath, filename)
        if os.path.exists(path):
            found = True
            break
        rootpath = os.path.split(rootpath)[0]
    if not found:
        raise MissingFileError(filename)
    return path


def get_file(label, folder):
    """Gets container file corresponding to container label in folder

    Parameters
    ----------
    label : str
        container label
    folder : str
        absolute path under which to look for container file

    Returns
    -------
    filename : str
        absolute path to file
    """
    accepted_extensions = ['', '.txt', '.tsv']
    for ext in accepted_extensions:
        fn = os.path.join(folder, label + ext)
        if os.path.exists(fn):
            break
    if not os.path.exists(fn):
        raise MissingFileError
    return fn


def get_array(fname, datatype, delimiter='\t'):
    """Returns Numpy structured array from text file

    Text file must be tab separated value and its columns must match the
    experiment descriptor file.

    Parameters
    ----------
    fname : str
        absolute path to text file to read
    datatype : Numpy readable datatype

    Returns
    -------
    numpy array
    """
    # big array of all cells
    arr = np.genfromtxt(fname, dtype=datatype, delimiter=delimiter)
    return arr


def datatype_parser(descriptor_file, sep=',', comment='!'):
    """Return Numpy datatype from descriptor file.

    Parameters
    ----------
    descriptor_file : str
        path to descriptor file
        should be a 2 columns file: col1: label, col2: type
    sep : str
        column separator
    comment : str
        comment character

    Returns
    -------
    datatype : list of couples ('label', type)
        datatype of flat files
    """
    datatypes = []
    abspath = os.path.abspath(os.path.expanduser(descriptor_file))
    with open(abspath, 'r') as f:
        for line in f.readlines():
            if line[0] != comment:
                seq = line.rstrip()
                if seq != '':
                    words = seq.split(sep)
                    datatypes.append(tuple(words[:2]))
    return datatypes


def has_container_data(abs_path):
    "Test whether absolute path has a sue"
    boo = False
    if 'containers' in os.listdir(abs_path):
        boo = True
    else:
        raise TextFileSystem(abs_path)
    return boo


def is_valid_experiment_folder(abs_path):
    res = glob.glob(abs_path)
    if not res:
        raise TextParsingError('No experiment folder at {}'.format(abs_path))
    res = glob.glob(os.path.join(abs_path, 'descriptor*'))
    if not res:
        raise TextParsingError('no descriptor file in experiment folder')
    res = glob.glob(os.path.join(abs_path, 'metadata*'))
    if not res:
        raise TextParsingError('no metadata file in experiment folder')
    return has_container_data(abs_path)


def find_containers(path):
    """Find container files

    Parameters
    ----------
    path : str, or pathlib.Path
        experiment root folder path

    Returns
    -------
    containers : list of pathlib.Path
    """
    if not isinstance(path, pathlib.Path):
        path = pathlib.Path(path)
    container_folder = path / 'containers'
    if not container_folder.exists():
        msg = 'Expected {}, not found. Check path.'.format(container_folder)
        raise TextParsingError(msg)
    containers = []
    for item in container_folder.iterdir():
        if item.is_file() and 'readme' not in item.stem.lower() and item.name[0] != '.':
            containers.append(item)
    return sorted(containers)


def find_metadata(path):
    """Find metadata associated to experiment pointed by path

    Parameters
    ----------
    path : str, or pathlib.Path
        experiment root folder path

    Returns
    -------
    meta : tunacell.io.metadata.Metadata instance
    """
    meta = metadata.load_metadata(path)
    return meta


def find_datatype(path):
    """Find datatype associated to text files

    Parameters
    ----------
    path : str, or pathlib.Path
        experiment root folder path

    Returns
    -------
    """
    descriptor_file = _check_up('descriptor.csv', path, 2)
    datatype = datatype_parser(descriptor_file)
    return datatype


class ContainerArrayParsingError(Exception):
    pass


class CellIdentifierError(ContainerArrayParsingError):
    """Class for Identifier Error"""
    pass


class CellParentError(ContainerArrayParsingError):
    """Class for parent identifier Error"""
    pass


def build_cells(arr, container=None, report_NaNs=True,
                extend_observables=False):
    """Read and store :class:`Cell` instances from structured text files).

    Text file must be tab separated value and its columns must match the
    experiment descriptor file.
    TODO: read descriptor from header?

    Parameters
    ----------
    arr : Numpy structured array
    container : :class:`Container` instance
    report_NaNs : boolean {True, False}
        whether to report for NaNs in text file
    extend_observables : boolean {False, True}
        whether to try to compute usual secondary observables such as age,
        volume, concentration...

    Returns
    -------
    list of :class:`Cell` instances
       Information is stored in attributes:
           * :attr:`bpointer`: backwards pointer, to parent cell
           * :attr:`data`: data as structured array
    """
    cells = []
    # big array of all cells
    if extend_observables:
        try:
            arr = compute_secondary_observables(arr)
        except ValueError as ve:
            msg = ('Extend observable failed, keep original array.\n'
                   '{}'.format(ve))
            logger.debug(msg)
    # store cellIDs
    cellIDs = np.unique(arr['cellID'])
    # when arr has got more than 1 frame
    if len(arr.shape) > 0:
        breaks = []  # where to split array
        previous_id = arr['cellID'][0]  # first cid
        for index, cid in enumerate(arr['cellID']):
            if cid != previous_id:
                previous_id = cid
                breaks.append(index)
        arrs = np.split(arr, breaks)
    # otherwise there's only one cell with a single frame
    else:
        arrs = [arr, ]

    del arr

    if report_NaNs:
        nan_labels = {}
    for arr in arrs:
        # first check that identifiers are unique
        cids = np.unique(arr['cellID'])
        if len(cids) > 1:
            raise CellIdentifierError('ids found: {}'.format(cids))
        else:
            cid = cids[0]
        pids = np.unique(arr['parentID'])
        if len(pids) > 1:
            msg = 'From cellID {}; found these parentIDs: {}'.format(cid, pids)
            raise CellParentError(msg)
        else:
            pid = pids[0]
        # create Cell instance and update bpointer when pid is valid
        cell = Cell(identifier=cid, container=container)
        if pid in cellIDs:  # register parent if and only if present as a recorded cell
            cell.bpointer = pid
        # record if NaN values appear
        if report_NaNs:
            for label, (dtype, offset) in arr.dtype.fields.items():
                # NaNs are implemented as np.nan for float types,
                if 'f' in dtype.kind:
                    if np.isnan(arr[label]).any():
                        if label not in nan_labels.keys():
                            nan_labels[label] = [cid, ]
                        else:
                            nan_labels[label].append(cid)
                # for integer types, they seem to be replaced by largest value
                elif ('u' in dtype.kind) or ('i' in dtype.kind):
                    if np.amax(arr[label]) == np.iinfo(dtype).max:
                        if label not in nan_labels.keys():
                            nan_labels[label] = [cid, ]
                        else:
                            nan_labels[label].append(cid)
        # attach data to Cell instance
        cell.data = arr
        cells.append(cell)
    if report_NaNs:
        msg = ('In container {}, found NaNs in following columns '.format(container.label) + ''
               '{}'.format(', '.join(nan_labels.keys())))
        logger.debug(msg)
    return cells
