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


def container_filename_parser(abs_path):
    "Parse experiment directory to retrieve FOV filenames."
    boo = is_valid_experiment_folder(abs_path)
    if boo:
        basenames = os.listdir(os.path.join(abs_path, 'containers'))
    # remove nondata files (so far remove the readme file)
    valids = []
    for basename in basenames:
        if 'readme' in basename or 'metadata' in basename:
            continue
        # skip hidden files
        if basename[0] == '.':
            continue
        valids.append(basename)
    return valids
#    filenames = [os.path.join(repo, bn) for bn in basenames]
#    return filenames


