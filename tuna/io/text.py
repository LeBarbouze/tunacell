#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Modules that defines import from/export to text files
"""
from __future__ import print_function

import os
import re
import glob

import numpy as np

from tuna.observable import Observable, FunctionalObservable


class TextParsingError(Exception):
    """General class for Errors while parsing text flat files."""
    pass


class MissingFileError(TextParsingError):
    """Class for missing file."""
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


class MissingFolderError(TextParsingError):

    def __init__(self, missing):
        self.missing = missing

    def __str__(self):
        msg = 'Missing folder: {}'.format(self.missing)
        return msg


# %% EXPORTING ANALYSIS FILES/FOLDERS

def get_analysis_path(exp, user_abspath=None, write=True):
    """Returns path to analysis folder.

    Parameters
    ----------
    exp : :class:`Experiment` instance
    user_abspath: str (default None)
        if given, will search within this path

    Returns
    -------
    analysis_path : str
        path to analysis folder corresponding to exp
    """
    # user defined main folder
    if user_abspath is not None:
        analysis = os.path.join(user_abspath, 'analysis')
    # canonical analysis folder
    else:
        # text
        if exp.filetype == 'text':
            analysis = os.path.join(exp.abspath, 'analysis')
        elif exp.filetype == 'h5':
            folder_up = os.path.split(exp.abspath)[0]
            analysis = os.path.join(folder_up, 'analysis')
        elif exp.filetype == 'simu':
            analysis = os.path.join(os.path.expanduser('~'), 'tmptuna',
                                    exp.label, 'analysis')
    if write and not os.path.exists(analysis):
        os.makedirs(analysis)
    return analysis


def _get_collections(path, basename='filterset'):
    """Build dict of folders with name starting with basename
    """
    if not os.path.exists(path):
        raise MissingFolderError('{} is missing'.format(path))
    p = re.compile(basename + '_(\d+)')
    ls = os.listdir(path)
    collec = {}  # dic repr(filt): (index, folder_path)
    # loop through directories and inspect each filterset directory
    for item in ls:
        path_to_item = os.path.join(path, item)
        if os.path.isdir(path_to_item):
            m = p.match(item)
            if m:
                sindex, = m.groups()
                index = int(sindex)
                rep = ''
                with open(os.path.join(path_to_item, item + '.txt'), 'r') as f:
                    rep = f.readline().rstrip()
                if rep:
                    collec[rep] = (index, path_to_item)
    return collec


def _get_item_path(folder, item, kind='filterset', write=True):
    """Returns path to corresponding item

    Parameters
    ----------
    folder : str
        absolute path of folder to look in
    item : :class:`FilterSet` to look for in folder
    kind : str {'filterset', 'condition'}
    write : bool {True, False}
        whether to write corresponding path on disk when it does not exist yet

    Returns
    -------
    index, path
    index : int
        integer label associated to item
    path : str
        absolute path to item
    """
    collec = _get_collections(folder, basename=kind)
    try:
        index, path = collec[repr(item)]
    except KeyError:
        used_indices = [c[0] for key, c in collec.items()]
        index = _get_new_index(used_indices, start_index=1)
        basename = '{}_{:02d}'.format(kind, index)
        if item.label is not None:
            basename += '_{}'.format(item.label)
        path = os.path.join(folder, basename)
        if write:
            os.makedirs(path)
            # write text file for filter description
            text_file = os.path.join(path, basename + '.txt')
            with open(text_file, 'w') as f:
                f.write('{}\n\n{}'.format(repr(item), str(item)))
    return index, path


def get_filters(analysis_path):
    return _get_collections(analysis_path, 'filterset')


def get_conditions(filter_path):
    return _get_collections(filter_path, 'condition')


def get_filter_path(analysis_path, fset, write=True):
    return _get_item_path(analysis_path, fset, kind='filterset', write=write)


def get_condition_path(filter_path, condition, write=True):
    # specific case for master : no further filter
    if condition is None or condition == 'master':
        path = os.path.join(filter_path, 'master')
        if write and not os.path.exists(path):
            os.makedirs(path)
        return 0, path
    else:
        return _get_item_path(filter_path, condition, kind='condition',
                              write=write)


def get_observable_path(condition_path, obs, write=True):
    if not os.path.exists(condition_path):
        raise MissingFolderError('condition-folder')
    basename = obs.name
    path = os.path.join(condition_path, basename)
    if write and not os.path.exists(path):
        os.makedirs(path)
        text_file = os.path.join(path, basename + '.txt')
        with open(text_file, 'w') as f:
            if isinstance(obs, Observable):
                f.write('{}\n\n{}\n\n{}\n\n{}'.format(repr(obs),
                                                      str(obs),
                                                      obs.as_latex_string(),
                                                      obs.as_string_table()))
            elif isinstance(obs, FunctionalObservable):
                f.write('{}'.format(basename))
    return path


def get_biobservable_path(condition_path, obss, write=True):
    if not os.path.exists(condition_path):
        raise MissingFolderError('condition-folder')
    basename = '---'.join([obs.name for obs in obss])
    path = os.path.join(condition_path, basename)
    if write and not os.path.exists(path):
        os.makedirs(path)
        # no writing of text dile description since univariate analysis did it
    return path

#def find_filterset_path(analysis_path, fset):
#    """Returns path corresponding to filterset.
#
#    Parameters
#    ----------
#    analysis_path : str
#        analysis folder
#    fset : FilterSet instance
#        filterset to be search within analysis folder
#
#    Returns
#    -------
#    (path, indices) when found
#
#    path : str if filter set is found, None if not found
#        absolute path to filterset folder
#    indices : list of int
#        list of used integers as indices for filtersets
#    """
#    path = None
#    p = re.compile('filterset_(\d+)')
#    if not os.path.exists(analysis_path):
#        raise MissingFolderError('analysis')
#    ls = os.listdir(analysis_path)
#    busy_indices = []
#    # loop through directories and inspect each filterset directory
#    for item in ls:
#        path_to_item = os.path.join(analysis_path, item)
#        if os.path.isdir(path_to_item):
#            m = p.match(item)
#            if m:
#                sindex, = m.groups()
#                index = int(sindex)
#                busy_indices.append(index)
#                ptf = os.path.join(path_to_item, item + '.txt')
#                with open(ptf, 'r') as f:
#                    line = f.readline()  # first line is repr(FilterSet)
#                    comp = line.rstrip()
#                    if comp == repr(fset):
#                        path = path_to_item
#    return (path, busy_indices)
#
#
#def find_observable_path(filterset_path, obs):
#    if not os.path.exists(filterset_path):
#        raise MissingFolderError('filterset')
#    path = os.path.join(filterset_path, obs.label)
#    return path
#
#
#def get_analysis_single(exp, fset, obs, user_abspath=None):
#    """Get analysis folder corresponding to fset, obs for exp.
#
#    Parameters
#    ----------
#    exp : Experiment instance
#    fset : FilterSet instance
#    obs : Observable instance
#    user_abspath : str (default None)
#        when user wants to override exp.abspath
#
#    Returns
#    -------
#    path to observable analysis folder
#
#    Raises
#    ------
#    MissingFolderError
#        with the missing level: 'filterset', 'observable'
#    """
#    res = None
#    analysis = get_analysis_path(exp, user_abspath=user_abspath)
#    path, busy_indices = find_filterset_path(analysis, fset)
#    if path is None:
#        raise MissingFolderError('filterset')
#    res = find_observable_path(path, obs)
#    if not os.path.exists(res):
#        raise MissingFolderError('observable')
#    return res
#
#
#def set_analysis_single(exp, fset, obs, user_abspath=None):
#    """Set up analysis folders for analysis of the dynamics.
#
#    Parameters
#    ----------
#    exp : Experiment instance
#    fset : FilterSet instance
#    obs : Observable instance
#    user_abspath : str (default None)
#        when user wants to override exp.abspath
#
#    Returns
#    -------
#    path to observable analysis folder
#
#    Notes
#    -----
#    The function is called when saving data is called upon
#    """
#    res = None
#    analysis = get_analysis_path(exp, user_abspath=user_abspath)
#    if not os.path.exists(analysis):
#        os.makedirs(analysis)
#    path, busy_indices = find_filterset_path(analysis, fset)
#    if path is None:
#        index = _get_new_index(busy_indices, limit=100)
#        item = 'filterset_{:02d}'.format(index)
#        path = os.path.join(analysis, item)
#        if not os.path.exists(path):
#            os.makedirs(path)
#        with open(os.path.join(path, item + '.txt'), 'w') as f:
#            f.write(repr(fset))  # first line : repr(fset)
#            f.write('\n\n')  # one blank line
#            f.write(str(fset))  # human readable description (str(fset))
#    res = find_observable_path(path, obs)
#    if not os.path.exists(res):
#        os.makedirs(res)
#    return res


def _get_new_index(busy_indices, start_index=0, limit=100):
    """Returns lowest integer lower than limit not in busy_indices list.

    Parameters
    ----------
    busy_indices : list of int
    limit : int (default 100)

    Returns
    -------
    int

    Examples
    --------

        >>> get_new_indew([0, 1, 2, 3], limit=100)
        4
        >>> get_new_indew([0, 2, 3, 4], limit=100)
        1

    """
    i = start_index
    while i < limit:
        if i not in busy_indices:
            break
        i += 1
    return i
