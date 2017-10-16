#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module defines the :class:`Colony` class that handle the tree-like
structure made by dividing cells.
"""
import treelib

import random

from numpy.random import randint
from tunacell.base.lineage import Lineage


class ColonyError(Exception):
    pass


class Colony(treelib.Tree):
    """Extension of treelib.Tree to add tunacell specific methods/attributes.

    Decomposition in independent lineages, to facilitate statistical analysis
    of dynamics, is relevant.

    Parameters
    ----------
    tree : Colony, or treelib.Tree instance
    deep : {False, True}
        option to deep copy input tree
    container : Container instance
        Container object to which current Colony belongs

    Note
    ----
    The initialization uses treelib.Tree initialization, adding two
    attributes .container, and .idseqs (for decomposition in lineages)

    See also
    --------
    The library treelib_ that implements the tree/node structures

    .. _treelib: https://github.com/caesar0301/treelib/blob/master/treelib/tree.py
    """

    def __init__(self, tree=None, deep=False, container=None):
        treelib.Tree.__init__(self, tree=tree, deep=deep)
        self.container = container
        self.idseqs = None
        return

    def add_cell_recursive(self, cell):
        """Function to add nodes recursively to a tree.

        Parameters
        ----------
        cell : Cell instance
           must have .parent and .childs attributes up-to-date
        """
        self.add_node(cell, parent=cell.bpointer)
        # print 'added %s to parent %s'%(cell.identifier, cell.bpointer)
        for ch in cell.childs:
            if ch.bpointer == cell.identifier:
                # we keep information about parents/childs
                # but link only if the backpointer still points to cell
                # which may be changed by prefiltering method
                self.add_cell_recursive(ch)
        return

    def decompose(self, independent=True):
        """Decompose tree onto lineages, i.e. sequences of cells

        Parameters
        ----------
        independent : bool (default True)
           with this option, returns list of sequences, where each cid appears
           only once. It is thus suited to perform statistical analysis.

        Returns
        -------
        idseqs : list of sequences of cell identifiers composing each lineage
        """
        if not independent:
            idseqs = self.paths_to_leaves()
        else:
            nids = self.expand_tree(mode=self.DEPTH, key=_randomise)
            idseqs = []
            seq = []
            for nid in nids:
                seq.append(nid)
                if self.get_node(nid).is_leaf():
                    idseqs.append(seq)
                    seq = []
#        if self.idseqs is None:
#            self.idseqs = independent_cell_lineages(self, only_ids=True)
#        # these may be used to build associated Lineage objects
        self.idseqs = idseqs
        return idseqs

    def iter_lineages(self, filt=None, size=None, shuffle=False):
        """Iterates through lineages.
        """
        if self.idseqs is None:
            idseqs = self.decompose()
        else:
            idseqs = self.idseqs[:]
        if shuffle:
            random.shuffle(idseqs)
        if filt is None:
            from tunacell.filters.lineages import FilterLineageAny
            filt = FilterLineageAny()
        count = 0
        for idseq in idseqs:
            if size is not None and count > size - 1:
                break
            lin = Lineage(self, idseq)
            if filt(lin):
                count += 1
                yield lin
        return


def _randomise(param):
    return random.uniform(0, 1)
