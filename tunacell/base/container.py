#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module defines how each container file (FoV) is handled
"""
from __future__ import print_function

import os
import re
import random
import numpy as np

import logging

from tunacell.io import text

from tunacell.base.cell import Cell, filiate_from_bpointer
from tunacell.base.datatools import compute_secondary_observables
from tunacell.base.colony import Colony, build_recursively_from_cells
from tunacell.base.observable import Observable, set_observable_list
from tunacell.filters.main import FilterTRUE
from tunacell.filters.trees import FilterTree


logger = logging.getLogger(__name__)


class ParsingContainerError(Exception):
    pass


class Container(object):
    """General class that stores reconstructed data per container.

    Read list of cells, associate microscopy data, make filiation,
    reconstruct trees. Note that only text files are fully supported so far.

    Parameters
    ----------
    filename : str
        container filename
        if exp is given, it must be an element of exp.files
    exp : Experiment instance
    path :str
        relative/absolute path to find the text file on disk
        [used only when exp is not provided] (default None)
    filetype : str {'text', 'simu'}
        [used only when exp is not provided] (default 'text')
    period : float
        time interval between two successive frames

    Attributes
    ----------
    exp : Experiment
        experiment to which container belongs
    abspath : str
        path to container file on disk
    filetype :str {'text'}
        type of file
    label : str
        label of container
    datatype : Numpy.dtype instance
        provides the datatype of raw data stored in each Cell instance .data
        attribute
    metadata : Metadata instance
        metadata of current container
    log : str
        logging message
    errs : sequence of Exceptions
        list of errors that occurred when reading data
    cells : list of Cell instances
        this list is created when reading data
    trees : list of Colony instance
        this list is created when reading data
    period : float
        periodicity of time-lapse acquisition (unit must match 'time' unit
        in data files)

    Methods
    -------
    get_cells()
        returns cells that are nodes of each container tree
    get_colony(cid)
        return colony to which belongs cell with identifier <cid>
    iter_colonies(filt=None, size=None, shuffle=False)
        iterates through all colonies of container
    postfilter(filt=None, verbose=False)
        to be used if new filters must be applied after reading data
        (at first reading prefilter method is applied)
    read_data(self, build=True, prefilt=None, extend_observables=True,
              report_NaNs=True)
        Container initialization does not read data in datafiles. Indeed it can
        be used just for parsing all container files without reading, and
        parsing is then faster. This method reads data, with the `build` option
        specifying that trees are (or aren't) build. The `prefilt` option
        is for applying a filter when reading data, i.e. rejecting cells that
        do not verify the filter assumption(s), before trees are constructed.
    """

    def __init__(self, label, exp):

        self.exp = exp  # experiment to which current Container belongs
        self.filetype = exp.filetype  # to be read or given

        # we will determine the following attributes in the __init__
        self.label = label  # to be read
        self.datatype = exp.datatype
        self.data = None  # Table read from file
#        try:
#            self.metadata = exp.metadata.loc[label]
#        except KeyError:
#            self.metadata = exp.metadata.loc[exp.label]
        self.metadata = exp.metadata.from_container(self.label)

        # these attributes are set to empty lists, will be loaded by .read_data
        self.cells = []
        self.trees = []

        # acquisition period
        self.period = self.metadata.loc['period']

        # cases against filetype
        if self.filetype == 'text':
            # identify absolute path to container file
            folder = os.path.join(self.exp.abspath, 'containers')
            # remove extension and store label
            self.label = label
#            print(label)
            if exp is not None and label not in exp.containers:
                msg = ''
                msg += 'This filename: {}'.format(label)
                msg += 'does not belong to exp: {}'.format(exp.abspath)
                raise ParsingContainerError(msg)
            self.abspath = text.get_file(label, folder)

            # Inherits datatype from Experiment datatype
            self.datatype = exp.datatype

        elif self.filetype == 'simu':
            pass

        return

    def read_data(self, build=True, prefilt=None, extend_observables=False,
                  report_NaNs=True):
        """Read the damn file

        Parameters
        ----------
        prefilt : Filter instance
            used to prefilter cells when reading data
        extend_observables : bool (default False)
            computes secondary obs
        report_NaNs : bool (default=True)
            report when NaNs are found in data
        """
        self.cells = []
        self.trees = []

        # TEXT FILETYPE
        if self.filetype == 'text':
            # Read cells from file
            arr = text.get_array(self.abspath, self.datatype, delimiter='\t')

        elif self.filetype == 'simu':
            pass
        
        self.data = arr

        if build:
            self.cells = build_cells(arr, container=self,
                                     extend_observables=extend_observables,
                                     report_NaNs=report_NaNs)

        
            self._build(prefilt=prefilt)

        return

    def _build(self, prefilt=None):
        """Builds colonies from list of cells read from files.

        Parameters
        ----------
        prefilt : :class:`tunacell.filters.cells.FilterCell` instance
            used to filter Cell instances at reading
        """
        self.make_filiation()
        if prefilt is not None:
            self.prefilter(filt=prefilt)
        self.make_trees()
        return

    def make_filiation(self):
        """Build filiation between cells.

        This method links cells objects.
        """
        if self.cells is not None:
            filiate_from_bpointer(self.cells)
#            for cell in self.cells:
#                childs = []
#                for cc in self.cells:
#                    if cc.bpointer == cell.identifier:
#                        childs.append(cc)
#                        cc.parent = cell
#                        cc.set_division_event()
#                cell.childs = childs
        return

    def prefilter(self, filt=None, verbose=False):
        """Filter at the cell level.

        Parameters
        ----------
        filt : Filter instance
           has to accept Cell instance as argument
           has to be callable
           returns True when cell is acccepted, False when not

        Notes
        -----
        This method acts of self.cells: it removes cells, and updates
        properties of parent/daughters cells when appropriate.

        Apply .prefilter only BEFORE building trees.

        See also
        --------
        The .postfilter method that filters cells at the tree level: when trees
        have been built after reading from file, they can be reconstructed
        upon filtering with this other method.
        """
        erased = []
        if verbose:
            msg = 'Prior to filter, we have {} cells.'.format(len(self.cells))
            print(msg)
        # check for supplementary observables to be computed
        raw_obs, func_obs = set_observable_list(filters=[filt, ])

        # compute suppl obs for all cells
        if raw_obs:
            for cell in self.cells:
                for obs in raw_obs:
                    cell.build(obs)
        for cell in self.cells:
            if filt is not None:
                if not filt(cell):
                    erased.append(cell)
                    # make Colony.add_cell_recursive non functional
                    cell.bpointer = None
                    if cell.parent:
                        cell.parent.childs.remove(cell)
                    # make daughter cells new roots
                    for ch in cell.childs:
                        ch.bpointer = None
        if verbose:
            msg = '{} cells do not pass filter.'.format(len(erased))
            print(msg)
        for cell in erased:
            self.cells.remove(cell)  # otherwise would be considered root
        # clean-up actions for computing extra obs
        # extra obs computation depends on tree decomposition
        # this will be done in lineage.get_timeseries()
        for cell in self.cells:
            for obs in raw_obs:
                del cell._sdata[obs.label]
        if verbose:
            msg = 'After filtering, we get {} cells.'.format(len(self.cells))
            print(msg)
#        self.metadata.filters.append(repr(boofunc))
        return

    def make_trees(self):
        """Build trees from list of cells."""
        self.trees = build_recursively_from_cells(self.cells, container=self)
#        self.trees = []
#        for cell in self.cells:
#            if cell.bpointer is None:  # test whether cell is root
#                tree = Colony(container=self)
#                tree.add_cell_recursive(cell)
#                self.trees.append(tree)
        return

    # TODO : postfiltering does not work with filter involving observables
    def postfilter(self, filt=None, verbose=False):
        """Rebuild trees after filtering on cells.

        Parameters
        ----------
        filt : Filter instance
           has to accept Cell instance as argument
           has to be callable
           returns True when cell is acccepted, False when not

        Notes
        -----
        This method acts on self.trees, the list of trees build from cells
        read in import file.

        Apply .postfilter AFTER having built trees. It updates the list of
        trees accordingly.

        See also
        --------
        the .prefilter method that acts on self.cells before building trees.
        
        .. warning:: does not compute supplementary observables in filters

        """
        new_trees = []
        if verbose:
            print('Starting postfiltering\n')
            print('Prior to filtering: {0} trees.'.format(len(self.trees)))

        for t in self.trees:
            outliers = []
            ns = t.all_nodes()
            if filt is not None:
                if filt.exonerate_root:
                    ns.remove(t.get_node(t.root))  # root has no data, so no filter
                for n in ns:
                    if not filt(n):
                        outliers.append(n)

            # sort outliers from leaves to root
            soutliers = sorted(outliers,
                               key=lambda x: t.level(x.identifier),
                               reverse=True)

            # with such ordering, no need for recursive action
            for n in soutliers:
                for ch in n.childs:
                    chid = ch.identifier
                    st = t.remove_subtree(chid)
                    st = Colony(tree=st, deep=False, container=self)
                    new_trees.append(st)
                p = n.parent
                if p is not None:
                    p.childs.remove(n)
                t.remove_node(n.identifier)

        # modify the list of trees
        self.trees += new_trees

        # update cells
        self.cells = [cell for tree in self.trees for cell in tree.all_nodes()]

        if verbose and filt is not None:
            msg = 'Post-filtering on cells: '
            msg += '{0} trees.'.format(len(self.trees))
            print(msg)
#        self.metadata.filters.append(repr(boofunc))
        return

    def get_cells(self):
        return [cell for tree in self.trees for cell in tree.all_nodes()]

    def get_colony(self, cid):
        """Retrieve colony to which belongs Cell instance with identifier 'cid'

        Parameters
        ----------
        cid : cell identifier (usually str)
        """
        found = False
        for colony in self.trees:
            if colony.contains(cid):
                found = True
                break
        if found:
            return colony
        else:
            msg = "There's no colony corresponding to {}".format(cid)
            msg += " in this container {}".format(self.label)
            raise ParsingContainerError(msg)

    def iter_colonies(self, filter_for_colonies='from_fset', size=None, shuffle=False):
        """Iterates through (already constructed colonies).

        Parameters
        ----------
        filt : filters.trees.FilterTree instance
           only colonies valid under filt will be yielded
        size : int
           maximal number of colonies before stopping iteration
        shuffle : bool
        """
        if filter_for_colonies is 'from_fset':
            colony_filter = self.exp.fset.colony_filter
        elif filter_for_colonies is None or filter_for_colonies is 'none':
            colony_filter = FilterTRUE()
        elif isinstance(filter_for_colonies, FilterTree):
            colony_filter = filter_for_colonies
        elif isinstance(filter_for_colonies, FilterTRUE):
            colony_filter = filter_for_colonies
        else:
            raise ValueError('"filter_for_colonies" parameter not recognized')
        trees = self.trees[:]
        if shuffle:
            random.shuffle(trees)
        count = 0
        for colony in trees:
            if size is not None and count > size - 1:
                break
            if colony_filter(colony):
                count += 1
                yield colony
        return

    def __repr__(self):
        msg = 'Container root: {}\n'.format(self.abspath)
        if self.cells is not None:
            msg += 'Cells: {}\n'.format(len(self.cells))
#            msg += '\n'.join([repr(c) for c in self.cells])
        else:
            msg += 'NO CELL\n'
        if self.trees is not None:
            msg += 'Trees: {}\n'.format(len(self.trees))
        else:
            msg += 'NO TREE\n'
        msg += repr(self.metadata)
        return msg

    def write_raw_text(self, path='.'):
        """Write data to text files.

        This function should be equivalent to copying text files...
        """
        cells = self.get_cells()
        arrays = []
        for cell in cells:
            arrays.append(cell.data)
        array = np.concatenate(arrays)
        fn = os.path.join(path, self.label + '.txt')
        fmt = []
        p = re.compile('(\w)(\d+)')
        for key, value in self.datatype:
            m = p.search(value)
            if m:
                kind, size = m.groups()
                # strings
                if kind == 'S':
                    add = '%{}c'.format(size)
                # integers
                elif kind in ['u', 'i']:
                    add = '%d'
                else:
                    add = '%.8e'
            else:
                add = '%.8e'
            fmt.append(add)
        np.savetxt(fn, array, fmt=fmt, delimiter='\t')
        return


# READING
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
            cid = str(cids[0])  # map to string (immutable)
        pids = np.unique(arr['parentID'])
        if len(pids) > 1:
            msg = 'From cellID {}; found these parentIDs: {}'.format(cid, pids)
            raise CellParentError(msg)
        else:
            pid = str(pids[0])  # map to string (immutable)
        # create Cell instance and update bpointer when pid is valid
        cell = Cell(identifier=cid, container=container)
        if pid != '0':  # this is the code for first recorded cells
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
