#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    tunacell.io.supersegger
    ^^^^^^^^^^^^^^^^^^^^^^^^

    module to parse supersegger data as input for tunacell processing
"""

from scipy.io import loadmat
import numpy as np

import sys
if sys.version_info[0] < 3:
    import pathlib2 as pathlib
else:
    import pathlib

from tunacell.base.cell import Cell


def find_containers(path):
    """Builds container list from experiment absolute path

    Parameters
    ----------
    path : str, or pathlib.Path

    Returns
    -------
    containers : list of pathlib.Path
        paths to containers
    """
    if not isinstance(path, pathlib.Path):
        path = pathlib.Path(path).expanduser().absolute()
    containers = [item for item in path.glob('xy*') if item.is_dir()]
    return containers


def load_container(path):
    """Load data from container folder path

    Parameters
    ----------
    path : str, or pathlib.Path
        path to the container folder

    Returns
    -------
    dict
        dictionary generated by loading the Matlab clist.mat file generated by
        Supersegger, see `SuperSegger Wiki Clist`_

    .. _SuperSegger Wiki Clist: https://github.com/wiggins-lab/SuperSegger/wiki/The-clist-data-file
    """
    if not isinstance(path, pathlib.Path):
        path = pathlib.Path(path).expanduser().absolute()
    clist = path / 'clist.mat'
    return loadmat(str(clist))


def read_header(mat, which='def3D'):
    """Reads header for 3D data

    Parameters
    ----------
    mat : dict
        io.loadmat of the clist.mat
    which : str {'def3D', 'def'}
        which data definition to use

    Returns
    -------
    header : list of str
        list of raw observables
    """
    dd = mat[which]
    header = [dd[0, k][0] for k in range(dd.shape[1])]
    return header


def read_ids(mat):
    """Builds dict cell ID to mother ID

    Parameters
    ----------
    mat : dict
        io.loadmat of the clist.mat

    Returns
    -------
    dict
        cell ID: parent cell ID relationship
    """
    header = read_header(mat, which='def')
    idx = header.index('Cell ID')
    pidx = header.index('Mother ID')
    data_ids = mat['data'][:, [idx, pidx]].astype('int')
    return {i: j for (i, j) in data_ids}


def build_cell_data(data, id2pid, header, time_array=None, period=None):
    """Build cell data structured array from Matlab data

    Parameters
    ----------
    data : (nobs, nframes) ndarray
    id2pid : dict
        dictionary cellID: parentID
    header : list of str
        header of axis 0 of ndarray
    time_array : (nframes, ) ndarray (default None)
        array of sampling times
    period : float (default None)
        when time_array is left None, period is used to map acquisition frame
        number to sampling time array (first acquisition time is set to 0.)

    Returns
    -------
    arr : Numpy structured array
    """
    # get observable index of cellID
    idx = header.index('Cell ID')
    # restrict to valid entries
    where, = np.where(np.logical_not(np.isnan(data[idx, :])))
    reduced = data[:, where]
    # shape of valid frames
    nobs, nframes = reduced.shape
    # build array of evaluation time
    if time_array is None:
        if period is None:
            raise ValueError('provide at least one defined argument')
        time = period * where
    else:
        time = time_array[where]
    # convert ids to array of int
    ids = reduced[idx].astype('int')
    if len(np.unique(ids)) > 1:
        raise ValueError('multiple ids for single cell')
    cid = np.unique(ids)[0]
    pid = id2pid[cid]
    pids = pid * np.ones_like(ids)
    # names for structured array
    names = 'cellID,parentID,time'
    formats='int,int,float'
    arrays = [ids, pids, time]
    for i in range(nobs):
        if i == idx:
            continue
        names += ',{}'.format(header[i])
        formats += ',float'
        arrays.append(reduced[i, :])
    # build record array
    arr = np.core.records.fromarrays(arrays, names=names, formats=formats)
    return arr


def build_cells(mat, container):
    """Build list of Cell instances

    Parameters
    ----------
    mat : dict
        dict obtained from loadmat on clist file
    container : tunacell.base.container.Container instance

    Returns
    -------
    list of :class:`Cell` instances
       Information is stored in attributes:
           * :attr:`bpointer`: backwards pointer, to parent cell
           * :attr:`data`: data as structured array
    """
    cells = []
    # builds dict of cell Id to parent ID
    dict_ids = read_ids(mat)
    # loop over cells in data3D
    header3D = read_header(mat, which='def3D')
    data3D = mat['data3D']
    ncells, nobs, nframes = data3D.shape
    for k in range(ncells):
        arr = build_cell_data(data3D[k, :, :], dict_ids, header3D,
                              time_array=None,
                              period=container.period)
        cid = arr['cellID'][0]
        pid = arr['parentID'][0]
        new = Cell(identifier=cid, container=container)
        new.data = arr
        if pid in dict_ids:  # points to parent if and only if itself a recorded cell
            new.bpointer = pid
        cells.append(new)
    return cells
