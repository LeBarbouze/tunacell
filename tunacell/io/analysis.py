#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  8 15:44:58 2018

@author: joachim
"""
import os
import re
import yaml
import logging
import warnings

from tunacell.io.text import (TextParsingError, MissingFileError,
                              MissingFolderError, MismatchFileError,
                              CorruptedFileError)

from tunacell.base.observable import Observable, FunctionalObservable

# import a bunch of filters to be able to load them using eval
from tunacell.filters.main import (FilterAND, FilterOR, FilterNOT, FilterTRUE,
                               FilterSet)
from tunacell.filters.cells import (FilterCellAny,
                                FilterCellIDparity,
                                FilterCompleteCycle,
                                FilterCycleFrames,
                                FilterCycleSpanIncluded,
                                FilterData,
                                FilterDaughters,
                                FilterHasParent,
                                FilterLengthIncrement,
                                FilterObservableBound,
                                FilterSymmetricDivision,
                                FilterTimeInCycle)
from tunacell.filters.trees import (FilterTreeAny,
                                FilterTreeDepth,
                                FilterTreeTimeIntersect)
from tunacell.filters.lineages import (FilterLineageAny,
                                   FilterLineageData,
                                   FilterLineageLength,
                                   FilterLineageTimeBound,
                                   FilterLineageTimeIntersect,
                                   FilterLineageTimeLength,
                                   FilterLineageWithCellProperty)
from tunacell.filters.containers import (FilterContainerAny,
                                     FilterContainerMetadataEquals)


logger = logging.getLogger(__name__)





# EXPORTING ANALYSIS FILES/FOLDERS

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
        if exp.filetype == 'text' or exp.filetype == 'supersegger':
            analysis = os.path.join(exp.abspath, 'analysis')
        elif exp.filetype == 'h5':
            folder_up = os.path.split(exp.abspath)[0]
            analysis = os.path.join(folder_up, 'analysis')
        elif exp.filetype == 'simu':
            analysis = os.path.join(os.path.expanduser('~'), 'tmptunacell',
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


def get_conditions(obs_path):
    return _get_collections(obs_path, 'condition')


def get_filter_path(analysis_path, fset, write=True):
    return _get_item_path(analysis_path, fset, kind='filterset', write=write)


def get_condition_path(obs_path, condition, write=True):
    # specific case for master : no further filter
    if condition is None or condition == 'master':
        path = os.path.join(obs_path, 'master')
        if write and not os.path.exists(path):
            os.makedirs(path)
        return 0, path
    else:
        return _get_item_path(obs_path, condition, kind='condition',
                              write=write)


def get_observable_path(filter_path, obs, write=True):
    if not os.path.exists(filter_path):
        raise MissingFolderError('filter-folder')
    basename = obs.name
    path = os.path.join(filter_path, basename)
    if write and not os.path.exists(path):
        logger.debug('Creating path {}'.format(path))
        os.makedirs(path)
        text_file = os.path.join(path, basename + '.txt')
        with open(text_file, 'w') as f:
            if isinstance(obs, Observable):
                f.write('{}\n\n{}\n\n{}\n\ncodestring: {}'.format(repr(obs),
                                                      obs.as_latex_string,
                                                      obs.as_string_table(),
                                                      str(obs)))
            elif isinstance(obs, FunctionalObservable):
                names = ', '.join([arg.name for arg in obs.observables])
                msg = '{}: FunctionalObservable({})'.format(basename, names)
                msg += '\n\n'
                for var_obs in obs.observables:
                    msg += '{}\n'.format(repr(var_obs))
                msg.rstrip()
                f.write(msg)
                # save serialized function
                source_file = os.path.join(path, basename + '_source.txt')
                with open(source_file, 'w') as sf:
                    sf.write('{}'.format(obs.source_f))
    # force write
    elif write and os.path.exists(path):
        logger.debug('May write existing file folders')
    # read mode: check that Observable representations match
    elif (not write) and os.path.exists(path):
        text_file = os.path.join(path, basename + '.txt')
        with open(text_file, 'r') as f:
            read_repr = f.readline().strip()
        if read_repr != repr(obs):
            if isinstance(obs, Observable):
                logger.debug('Obs path {} does not match argument {}'.format(path, obs))
                logger.debug('Changing name by appending a letter')
                raise MismatchFileError(level='observable')
            elif isinstance(obs, FunctionalObservable):
                msg = ('Impossible to check whether FunctionalObservable '
                       'matches since the export of its combining function is '
                       'not set up.')
                warnings.warn(msg)
        else:
            logger.debug('Reading matching observable path {}'.format(path))
    return path


def get_biobservable_path(filter_path, obss, write=True):
    """Get folder for bivariate analysis

    Parameters
    ----------
    filter_path : str
        parent folder, should be a filterset path
    obss : couple of :class:`Observable` instances
    write : bool {True, False}
        whether to write new path or not

    Returns
    -------
    path : str
        path to the bivariate analysis folder
    """
    if not os.path.exists(filter_path):
        raise MissingFolderError('filter-folder')
    basename = '---'.join([obs.name for obs in obss])
    path = os.path.join(filter_path, basename)
    if write and not os.path.exists(path):
        os.makedirs(path)
        # no writing of text dile description since univariate analysis did it
    return path


def read_count_file(filter_path):
    """Read yaml file for count

    Parameters
    ----------
    filter_path : str
        path to a FilterSet folder

    Returns
    -------
    counts : dict
        keys are cells, lineages, colonies, containers
    """
    count_file = os.path.join(filter_path, '.counts.yml')
    if not os.path.exists(count_file):
        raise MissingFileError
    with open(count_file, 'r') as f:
        counts = yaml.load(f)
    # check that information is correctly stored
    a = 'cells' in counts
    b = 'lineages' in counts
    c = 'colonies' in counts
    d = 'containers' in counts
    if not (a and b and c and d):
        raise CorruptedFileError
    return counts


def write_count_file(filter_path, counts):
    """Write yaml file"""
    count_file = os.path.join(filter_path, '.counts.yml')
    with open(count_file, 'w') as f:
        yaml.dump(counts, stream=f, default_flow_style=False)


# PRINTING STUFF FROM ANALYSIS TEXT FILES

def _print_collections(parent_folder, kind='filterset'):
    """Print list of filtersets/conditions

    Parameters
    ----------
    parent_folder : str
        parent folder in which to look for filtersets/conditions
    kind : str {'filterset', 'condition'}
    """
    msg = 'Looking for {}s under {} ...'.format(kind, parent_folder)
    collec = _get_collections(parent_folder, basename=kind)
    # order items using index
    as_list = sorted([(index, path_to_item, rep)
                      for rep, (index, path_to_item) in collec.items()],
                     key=lambda x: x[0])
    if len(as_list) == 0:
        msg += '\n\n Nothing here. Move along.'
    for (index, path_to_item, rep) in as_list:
        basename = os.path.basename(path_to_item)
        consensus = '{}_{:02d}_(\S*)'.format(kind, index)
        chain = re.compile(consensus)
        m = chain.match(basename)
        name = ''
        if m:
            name, = m.groups()
        if not name:
            name = '(none)'
        fname = os.path.join(path_to_item, basename + '.txt')
        if not os.path.exists(fname):
            raise MissingFileError('Missing description file under {}'.format(path_to_item))
        rep, human = _read_first_remaining(fname)
        msg += '\n\n{}. name: {} path: {}'.format(index, name, path_to_item)
        msg += '\nrep: {}'.format(rep)
        if human:
            msg += '\n{}'.format(human)
    print(msg)
    print()
    return


def print_filtersets(exp):
    """Print out filtersets saved in analysis folder

    Parameters
    ----------
    exp : :class:`Experiment` instance
    """
    analysis_path = get_analysis_path(exp, write=False)
    if not os.path.exists(analysis_path):
        print('There is no analysis folder. Compute, export, and come back later')
    _print_collections(analysis_path, kind='filterset')
    return


def print_conditions(exp, fset, obs):
    """Print out conditions used for input observable

    Parameters
    ----------
    exp : :class:`Experiment` instance
    fset : :class:`FilterSet` instance
    obs : :class:`Observable` instance
    """
    analysis_path = get_analysis_path(exp, write=False)
    _, filter_path = get_filter_path(analysis_path, fset, write=False)
    obs_path = get_observable_path(filter_path, obs, write=False)
    _print_collections(obs_path, kind='condition')
    return


def print_observables(exp, fset):
    """Print out observables that have been analyzed

    Parameters
    ----------
    exp : :class:`Experiment` instance
    fset : :class:`FilterSet` instance
    """
    analysis_path = get_analysis_path(exp, write=False)
    _, filter_path = get_filter_path(analysis_path, fset, write=False)
    msg = 'Looking for observables under {} ...'.format(filter_path)
    items = os.listdir(filter_path)
    candidates = [item for item in items
                  if os.path.isdir(os.path.join(filter_path, item))]
    valids = []
    for name in candidates:
        abs_path = os.path.join(filter_path, name)
        fname = os.path.join(abs_path, name + '.txt')
        if os.path.exists(fname):
            rep, human = _read_first_remaining(fname)
            if 'Observable' in rep or 'FunctionalObservable' in rep:
                valids.append(name)
                msg += '\n\n{} path: {}'.format(name, abs_path)
                msg += '\nrep: {}'.format(rep)
                if human:
                    msg += '\n\n{}'.format(human)
    if len(valids) == 0:
        msg += 'Nothing there. Move along'
    print(msg)
    print()
    return


# LOADING STUFF FROM TEXT FILES

class ImpossibleToLoad(ValueError):
    pass


def load_item_from_path(path):
    """Returns an evaluated object from path"""
    basename = os.path.basename(path)
    fname = os.path.join(path, basename + '.txt')
    if not os.path.exists(fname):
        raise MissingFileError('Missing description file under {}'.format(path))
    rep, human_string = _read_first_remaining(fname)
    # so far, only FunctionalObservable are not loadable (TO FIX)
    if 'FunctionalObservable' in rep:
        raise ImpossibleToLoad('FunctionalObservable are not loadable')
    return eval(rep)


# other functions


def _read_first_remaining(filename):
    """Get first line of file, and remaining content as couple.

    Parameters
    ----------
    filename : str
        absolute path to file

    Returns
    -------
    (first line, remaining content): str, str
    """
    with open(filename, 'r') as f:
        first = f.readline()
        msg = ''
        for line in f.readlines():
            msg += line
    return (first.rstrip(), msg.lstrip().rstrip())


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

        >>> _get_new_index([0, 1, 2, 3], limit=100)
        4
        >>> _get_new_index([0, 2, 3, 4], limit=100)
        1

    """
    i = start_index
    while i < limit:
        if i not in busy_indices:
            break
        i += 1
    return i
