#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module defines the :class:`Colony` class that handle the tree-like
structure made by dividing cells.
"""
import treelib

import random

from numpy.random import randint
from tuna.base.lineage import Lineage


class ColonyError(Exception):
    pass


class Colony(treelib.Tree):
    """Extension of treelib.Tree to add tuna specific methods/attributes.

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
            from tuna.filters.lineages import FilterLineageAny
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


# %% CODE BELOW IS DEPRECATED

def find_super_leaf(tree):
    """Find leaf with maximal depth.

    Arguments
    ---------
    tree -- treelib.Tree instance

    Note
    ----
    super leaf is chosen randomly from the set of maximal depth leaves.
    """
    super_leaf_ids = []
    depth = tree.depth()
    if depth == 0:
        return tree.root
    super_leaf = None
    for l in tree.leaves():
        if tree.level(l.identifier) == depth:
            super_leaf_ids.append(l.identifier)
#            last_time = l.data[-1]['time']
#            if last_time > max_time:
#                super_leaf = l.identifier
#                max_time = last_time
    super_leaf = super_leaf_ids[randint(0, len(super_leaf_ids))]
    return super_leaf


# 2. Recursive function
def go_up(subtree, lines):
    """Function that build recursively independent lineages from a tree.

    Arguments
    ---------
    subtree -- treelib.Tree instance
    lines -- list of sequence of Cell identifiers
    """
    leaf = find_super_leaf(subtree)
    main_line = [nid for nid in subtree.rsearch(leaf)]
    lines.append(main_line)
    index = 0
    nid = main_line[index]
    while nid != subtree.root:
        sibs = subtree.siblings(nid)
        if sibs:
            s = sibs[0]  # at most one sibling
            go_up(treelib.Tree(subtree.subtree(s.identifier), deep=True), lines)
        index += 1
        subtree.remove_node(nid)
        nid = main_line[index]
    return


def independent_cell_lineages(tree,
                              remove_nodata=True,
                              remove_leaves=False,
                              only_ids=True):
    """Returns a list of independent lineages of tree.

    Argument
    --------
    tree -- treelib.Tree instance

    Parameters
    ----------
    remove_nodata -- boolean, default True
       remove cell id from output sequence when no data is associated to cell
    remove_leaves -- boolean, default False
       remove the last cell of output sequence
    only_ids -- boolean, default True
       if True, returns sequences of cell identifiers;
       if False, returns sequences of cell instances.

    Returns
    -------
    Independent sequences of cells from a tree structure.

    Notes
    -----
    Each lineage is sorted from root to leave.
    """
    idlines = []
    test_tree = treelib.Tree(tree, deep=True)
    go_up(test_tree, idlines)
    for line in idlines:
        if remove_leaves:
            line.remove(line[0])
        if remove_nodata:
            for cid in line:
                if tree.get_node(cid).data is None:
                    line.remove(cid)
    idlines_clean = [line[::-1] for line in idlines if line]
    if only_ids:
        # return list of cell labels (default)
        return idlines_clean
    else:
        # return list of Cell instances
        cell_lines = []
        for line in idlines_clean:
            cline = [tree.get_node(cid) for cid in line]
            cell_lines.append(cline)
        return cell_lines
