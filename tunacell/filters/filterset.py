#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""FilterSet defines filters for each type: cell/lineage/tree/container"""
from __future__ import print_function

try:
    from StringIO import StringIO  # python2
except ImportError:
    from io import StringIO  # python3

from tunacell.filters.cells import FilterCellAny
from tunacell.filters.lineages import FilterLineageAny
from tunacell.filters.trees import FilterTreeAny
from tunacell.filters.containers import FilterContainerAny


class FilterSet(object):
    """Collects filters of each type in a single object"""

    def __init__(
        self,
        label=None,
        filtercell=FilterCellAny(),
        filterlineage=FilterLineageAny(),
        filtertree=FilterTreeAny(),
        filtercontainer=FilterContainerAny(),
    ):
        """Instantiate a filter set.

        Parameters
        ----------
        label : str (default None)
            give a label to current FilterSet
        filtercell : :class:`FilterGeneral` instance
            must be of type 'CELL'
        filterlineage : :class:`FilterGeneral` instance
            must be of type 'LINEAGE'
        filtertree : :class:`FilterGeneral` instance
            must be of type 'TREE'
        filtercontainer : :class:`FilterGeneral` instance
            must be of type 'CONTAINER'
        """
        self.label = label
        # initialize filters
        self._obs = []  # needed to compute testable observables
        self.cell_filter = filtercell
        self.lineage_filter = filterlineage
        self.colony_filter = filtertree
        self.container_filter = filtercontainer
        self._collect_hidden_obs()

        self._named_list = [
            ("label", self.label),
            ("filtercell", self.cell_filter),
            ("filterlineage", self.lineage_filter),
            ("filtertree", self.colony_filter),
            ("filtercontainer", self.container_filter),
        ]
        return

    def _collect_hidden_obs(self):
        """Returns the list of hidden observables (needed for computations)"""
        self._obs = []
        for filt in [
            self.cell_filter,
            self.lineage_filter,
            self.colony_filter,
            self.container_filter,
        ]:
            for suppl_obs in filt._obs:
                if suppl_obs not in self._obs:
                    self._obs.append(suppl_obs)
        return

    @property
    def obs(self):
        """Provides the list of hidden observables"""
        return self._obs

    def __repr__(self):
        name = type(self).__name__
        chain = name + "("
        for (name, item) in self._named_list:
            chain += "{}={}, ".format(name, repr(item))
        chain += ")"
        return chain

    def __str__(self):
        label = ""
        if self.label is not None:
            label += "label: {}\n".format(self.label)
        for filt in [
            self.cell_filter,
            self.lineage_filter,
            self.colony_filter,
            self.container_filter,
        ]:
            if filt._type == "ANY":
                continue
            for index, line in enumerate(StringIO(str(filt))):
                if index == 0:
                    label += "* "
                else:
                    label += "  "
                label += line
            label += "\n"
        return label.rstrip()
