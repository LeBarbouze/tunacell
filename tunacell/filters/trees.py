#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
tunacell package
============

filters/trees.py module
~~~~~~~~~~~~~~~~~~~~~~~~~~

Classes to filter trees.
"""
import numpy as np
import warnings

from tunacell.filters.main import FilterGeneral, bounded, intersect


class FilterTree(FilterGeneral):
    "General class for filtering tree objects (treelib.Tree instances)"

    _type = 'TREE'


class FilterTreeAny(FilterTree):

    def __init__(self):
        self.label = 'Always True'
        return

    def func(self, tree):
        return True


class FilterTreeDepth(FilterTree):

    def __init__(self, lower_bound=3, upper_bound=None):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        label = 'Tree depth filter: '
        label += '{0} <= tree_depth <= {1}'.format(lower_bound, upper_bound)
        self.label = label
        return

    def func(self, tree):
        import treelib
        boo = False
        if isinstance(tree, treelib.Tree):
            boo = bounded(tree.depth(),
                          lower_bound=self.lower_bound,
                          upper_bound=self.upper_bound)
        else:
            warnings.warn('Argument is not a tree...')
        return boo


class FilterTreeTimeIntersect(FilterTree):

    def __init__(self, lower_bound=None, upper_bound=None):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound
        label = '{} < tree time span < {}'.format(lower_bound, upper_bound)
        self.label = label
        return

    def func(self, tree):
        root = tree.get_node(tree.root)
        values = [np.amin(root.data['time']), ]
        tmaxs = []
        for leaf in tree.leaves():
            if leaf.data is not None:
                tmaxs.append(np.amax(leaf.data['time']))
        if tmaxs:
            values.append(np.amax(tmaxs))
        return intersect(values, lower_bound=self.lower_bound,
                         upper_bound=self.upper_bound)
