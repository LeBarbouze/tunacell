#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module is designed to inspect experiment folders to look for stored
analysis files.
"""
from __future__ import print_function
import os
import re
import warnings

import treelib

from tuna.io.text import is_valid_experiment_folder
from tuna import Observable
from tuna.parser import Parser

from tuna.filters.main import (FilterAND, FilterOR, FilterNOT, FilterTRUE,
                               FilterSet)
from tuna.filters.cells import (FilterCellAny,
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
from tuna.filters.trees import (FilterTreeAny,
                                FilterTreeDepth,
                                FilterTreeTimeIntersect)
from tuna.filters.lineages import (FilterLineageAny,
                                   FilterLineageData,
                                   FilterLineageLength,
                                   FilterLineageTimeBound,
                                   FilterLineageTimeIntersect,
                                   FilterLineageTimeLength,
                                   FilterLineageWithCellProperty)
from tuna.filters.containers import (FilterContainerAny,
                                     FilterContainerMetadataEquals)
try:
    from StringIO import StringIO  # python2
except ImportError:
    from io import StringIO  # python3


# %% ANALYSIS FOLDER SNIFFING
class SniffingError(Exception):
    pass


(EXPERIMENT, ANALYSIS, FILTERSET, OBSERVABLE, CONDITION) = list(range(5))


class Sniffer(object):
    """A class to explore the experiment folder"""

    def __init__(self, path=None):
        self.path = None
        if path is not None:
            self.path = os.path.abspath(os.path.expanduser(path))
            if not is_valid_experiment_folder(self.path):
                msg = ('Not a valid Experiment folder.\n'
                       'Use .sniff(path=<your path>)')
                warnings.warn(msg)
        if self.path is None:
            warnings.warn('Use .sniff(<your_path>) to load and sniff an experiment.')
        else:
            _, self.label = os.path.split(self.path)
        self._tree = None
        self.sniff()
        return

    def sniff(self, path=None):
        if self.path is None:
            if path is None:
                msg = ('Load an experiment folder path as argument please.')
                warnings.warn(msg)
                return
            else:
                self.path = path
#        analysis = os.path.join(self.path, 'analysis')
#        if not os.path.exists(analysis):
#            msg = ('No analysis folder under {}'.format(self.path))
#            warnings.warn(msg)
#            return
        fold = SniffFolder(self.path)
        root = treelib.Node(tag=self.label, data=fold)
        root.level = 0  # necessary for later purposes
        tree = treelib.Tree()
        tree.add_node(root, parent=None)
        _add_children_nodes(tree, root)
        self._tree = tree
        tree.show()
        return

    def show(self):
        if self._tree is not None:
            self._tree.show()
        else:
            print('No build tree')
        return

    def show_details(self, from_node=None, min_level=0, max_level=None):
        """Print text output with description of each node in tree

        Parameters
        ----------
        from_node : :class:`Node` instance
            folder tree is parsed from this node.
            :attr:`data` must be a SniffFolder instance
        min_level : int (default 0)
            minimum level to look up in folder tree, used only when
            :param:`from_node` is left to None
        max_level : int (default None)
            maximum level to look up in folder tree
        """
        if from_node is not None:
            starting_nodes = [from_node, ]
        else:
            starting_nodes = []
            for node in self._tree.all_nodes():
                if not hasattr(node, 'level'):
                    level = self._tree.level(node.identifier)
                else:
                    level = node.level
                if level == min_level:
                    starting_nodes.append(node)
        msg = ''
        if max_level is not None:
            max_depth = max_level
        else:
            max_depth = self._tree.depth()
        for starting_node in starting_nodes:
            for nid in self._tree.expand_tree(nid=starting_node.identifier):
                node = self._tree.get_node(nid)
                depth = self._tree.level(node.identifier)
                if depth > max_depth:
                    continue
                # absolute path for experiment, relative path otherwise
                if depth == 0:
                    path = node.data.path
                else:
                    _, last = node.data.path.split(self.label)
                    path = self.label + last
                indent = '  ' * depth
                loc = ('{}label: {}'.format(indent, node.data.label) + '\n'
                       '{}path: {}'.format(indent, path) + '\n'
                       '{}description:'.format(indent) + '\n')
                for line in StringIO(node.data.human).readlines():
                    loc += '{}{}'.format(indent, line)
                loc += '\n\n'
                msg += loc
        if not starting_nodes:
            msg = ('No node at level {}'.format(min_level))
        print(msg.rstrip())
        return

    def show_filtersets(self, from_node=None):
        self.show_details(from_node=from_node,
                          min_level=FILTERSET,
                          max_level=FILTERSET)
        return

    def show_obs(self, from_node=None):
        self.show_details(from_node=from_node,
                          min_level=OBSERVABLE,
                          max_level=OBSERVABLE)
        return

    def show_conditions(self, from_node=None):
        self.show_details(from_node=from_node,
                          min_level=CONDITION,
                          max_level=CONDITION)
        return

    def get_main_node(self, label='analysis'):
        """Returns node corresponding to experiment main folder name label."""
        res = None
        root = self._tree.get_node(self._tree.root)
        for nid in root.fpointer:
            node = self._tree.get_node(nid)
            if node.data.label == label:
                res = node
        if res is None:
            raise SniffingError('No folder label {}'.format(label))
        return res

    def get_filterset_node(self, index, folder='analysis'):
        """Returns filterset corresponding to index."""
        res = None
        level = FILTERSET
        from_node = self.get_main_node(label=folder)
        for nid in from_node.fpointer:
            node = self._tree.get_node(nid)
            if self._tree.level(nid) == level and node.data.index == index:
                res = node
        if res is None:
            raise SniffingError('No filterset_{:02d}'.format(index))
        return res

    def get_obs_node(self, codestring, filterset_index, folder='analysis'):
        """Returns obs node corresponding to observable codestring.

        Parameters
        ----------
        codestring : str
            codestring encoding :class:`Observable` parameters
        filterset_index: int
            sets the filterset for which observable has been computed
        folder : str
            this is the main folder under which filterset are defined

        Returns
        -------
        :class:`Node` instance
        """
        res = None
        fnode = self.get_filterset_node(filterset_index, folder=folder)
        for nid in fnode.fpointer:
            node = self._tree.get_node(nid)
            if node.data.label == codestring:
                res = node
        if res is None:
            raise SniffingError('No obs codestring {}'.format(codestring))
        return res

    def get_condition_node(self, condition_index, codestring, filterset_index,
                           folder='analysis'):
        """Returns condition node corresponding to observable codestring

        Parameters
        ----------
        condition_index : int
            condition index to get. Note that 0 corresponds to 'master';
            nodes corresponding to condition Filters start at index 1
        codestring : str
            observable codestring
        filterset_index : int
            filterset index
        folder : str (default 'analysis')
            main folder to look into

        Returns
        -------
        :class:`Node` instance
            node instance of folder tree, with :attr:`.data` pointing
            to corresponding :class:`SniffingFolder` instance.
        """
        onode = self.get_obs_node(codestring, filterset_index, folder)
        res = None
        for cid in onode.fpointer:
            node = self._tree.get_node(cid)
            if condition_index == 0 and node.data.label == 'master':
                res = node
                break
            elif node.data.index == condition_index:
                res = node
                break
        if res is None:
            raise SniffingError('No condition_{:02d}'.format(condition_index))
        return res

    def get_node(self, path):
        """Returns node corresponding to SniffFolder instance at path=
        """
        res = None
        for node in self._tree.all_nodes():
            if _filter_path(path, node):
                res = node
        if res is None:
            raise SniffingError('No subfolder {}'.format(path))
        return res

    def load_filterset(self, filterset_index=0, folder='analysis',
                       abspath=None):
        """Return chosen :class:`FilterSet`

        Two modes of search: either using absolute path, either using the
        filterset index and the main folder.

        Parameters
        ----------
        filterset_index : int (default 0)
            index of the filterset to load
        folder : str (default 'analysis')
            main folder under which to look for filterset
        abspath : str (default None)
            when given, try to open the filterset using the absolute path,
            otherwise uses the other keyword arguments

        Returns
        -------
        :class:`FilterSet` instance

        Raises
        ------
        :exception:`SnifferError`
            when searching for given parameters fails
        """
        res = None
        if abspath is not None:
            node = self.get_node(abspath)
            if node.data.level != self.FILTERSET:
                raise SniffingError('This is no filterset')
            res = eval(node.data.representation)
        else:
            node = self.get_filterset_node(filterset_index, folder=folder)
            res = eval(node.data.representation)
        return res

    def load_obs(self, codestring='codestring',
                 filterset_index=0, folder='analysis',
                 abspath=None):
        """Returns chosen :class:`Observable` instance.

        Selection via absolute path, or via other keyword arguments.

        Parameters
        ----------
        codestring : str (default 'codestring')
            this is the observable codestring representation to load.
            Use default only when abspath is given.
        filterset_index : int (default 0)
            index of the filterset to be loaded.
        folder  : str (default 'analysis')
            main folder under which filterset and observable are looked upon,
            no incidence when abspath is given.
        abspath : str (default None)
            when None, use the other keyword arguments to load observable;
            otherwise use only the given absolute path

        Returns
        -------
        :class:`Observable` instance

        Raises
        ------
        :exception:`SniffingError`
            when it is not possible to find the required observable, a small
            informative message should indicate which part of the search went
            wrong.
        """
        res = None
        if abspath is not None:
            node = self.get_node(abspath)
            if node.data.level != OBSERVABLE:
                raise SniffingError('This is no observable')
            res = eval(node.data.representation)
        else:
            node = self.get_obs_node(codestring, filterset_index,
                                     folder=folder)
            res = eval(node.data.representation)
        return res

    def load_condition_element(self, condition_index=1,
                               codestring='codestring',
                               filterset_index=0,
                               folder='analysis',
                               abspath=None):
        """Returns the :class:`FilterGeneral` instance defining condition.

        Parameters
        ----------
        condition_index : int (default 1)
            index of the condition to be loaded. Avoid 0 as it is reserved for
            'master' (unconditioned data)
        codestring : str (default 'codestring')
            observable codestring
        filterset_index : int (default 0)
            index of the chosen filterset
        folder : str (default 'analysis')
            main folder under which to look for

        Returns
        -------
        :class:`FilterGeneral` instance
            filter that defines the condition

        Raises
        ------
        :exception:`SniffingError`
            when search is unsuccessful
        """
        res = None
        if abspath is not None:
            node = self.get_node(abspath)
            if node.data.level != CONDITION:
                raise SniffingError('This is no condition')
            res = eval(node.data.representation)
        else:
            node = self.get_condition_node(condition_index, codestring,
                                           filterset_index, folder)
            res = eval(node.data.representation)
        return res

    def load_conditions(self, codestring='codestring', filterset_index=0,
                        folder='analysis',
                        abspath=None):
        """Returns list of :class:`FilterSet` instances defined under obs node.
        
        They are the set of conditions user applied on saved files.

        Parameters
        ----------
        codestring : str (default 'codestring')
            leave to default only when abspath is given
        filterset_index : int (default 0)
            index of the chosen filterset
        folder : str (default 'analysis')
            main folder under which to look for
        abspath : str (default None)
            if left None, other keyword arguments are used

        Returns
        -------
        list of :class:`FilterSet` instances
        """
        node = None
        if abspath is not None:
            node = self.get_node(abspath)
            if node.data.level != OBSERVABLE:
                msg = '{} should point to observable'.format(abspath)
                raise SniffingError(msg)
        else:
            node = self.get_obs_node(codestring, filterset_index, folder)
        if node is None:
            raise SniffingError('This is not a node')
        selections = []
        for cid in node.fpointer:
            cnode = self._tree.get_node(cid)
            if cnode.data.label != 'master':
                selections.append(eval(cnode.data.representation))
        return selections


def load_framework(sniffer, filterset_index=0, codestring=None,
                   folder='analysis'):
    """Loads a saved framework.

    This function sets :class:`Parser` (as a combination of both
    :class:`Experiment` and :class:`Filterset`), and :class:`Observable`
    and associated computed list of :class:`Filterset` instances when
    observable codestring is given.

    Parameters
    ----------
    filterset_index : int (default 0)
        index of filterset to be loaded
    codestring : str (default None)
        Observable codestring to be loaded
    folder : str (default 'analysis')
        folder under which filterset_index is search for

    Returns
    -------
    (parser, observable, conditions)
    parser : :class:`Parser` instance
        defines how to parse data, including experiment and filter set
    observable : :class:`Observable` instance
        observable that was analyzed
    conditions : list of :class:`FilterSet` instances
        the various conditions defined in files

    See also
    --------
    :class:`Sniffer`

    Warning
    -------
    only text experiment files can be loaded at this stage
    """
    parser = Parser()
    parser.load_experiment(sniffer.path, filetype='text')
    fnode = sniffer.get_filterset_node(filterset_index, folder=folder)
    fset = eval(fnode.data.representation)
    parser.fset = fset
    obs, cset = None, None
    if codestring is not None:
        onode = sniffer.get_obs_node(codestring, filterset_index,
                                     folder=folder)
        obs = eval(onode.data.representation)
        selections = []
        for nid in onode.fpointer:
            node = sniffer._tree.get_node(nid)
            rep = node.data.representation
            if node.data.label != 'master' and rep is not None:
                selections.append(eval(rep))
        cset = selections
    return parser, obs, cset


class SniffFolder(object):
    """Excludes experiment root file that is parsed in Sniffer"""

    def __init__(self, path):
        self.path = path
        parent, label = os.path.split(path)
        self.label = label
        self.level = None  # exp, analysis, filterset, obs, condition
        self.representation = None
        self.human = None
        self._contains = []  # sub-folders to explore
        self._sniff()  # will check what kind of folder this is representing
        return

    def _sniff(self):
        """Snif according to label.
        """
        fpattern = re.compile('filterset_(\d+)')
        opattern = re.compile('T([a-z])(\d*[\.,]*\d*)M([a-z\-]+)J(\d+)')
        cpattern = re.compile('condition_(\d+)')
        # parse files under path, get only directories
        self._contains = [item for item in os.listdir(self.path)
                          if os.path.isdir(os.path.join(self.path, item))]
        fname = None
        self.level = None  # default
        if self.label in ['analysis', 'master']:
            if self.label == 'analysis':
                self.level = ANALYSIS
            elif self.label == 'master':
                self.level = CONDITION
            return
        elif fpattern.match(self.label):
            self.level = FILTERSET
            sindex, = fpattern.match(self.label).groups()
            self.index = int(sindex)
            fname = os.path.join(self.path, self.label + '.txt')
        elif opattern.search(self.label):
            self.level = OBSERVABLE
            self.codestring = self.label
            fname = os.path.join(self.path, 'observable.txt')
        elif cpattern.match(self.label):
            self.level = CONDITION
            sindex, = cpattern.match(self.label).groups()
            self.index = int(sindex)
            fname = os.path.join(self.path, self.label + '.txt')

        if fname is not None:
            self.representation, self.human = _read_first_remaining(fname)
        return


def _add_children_nodes(tree, parent):
    """Iteratively add filterset sniffers

    Parameters
    ----------
    tree : :class:`treelib.Tree` instance
        tree that recapitulates folder organization and sniffing content
    parent : :class:`treelib.Node` instance
        parent node (usually analysis node) to which filterset nodes will be
        associated

    Yields
    ------
    fnode : :class:`treelib.Node` instance
        child node to parent, containing
    """
    for item in parent.data._contains:
        path = os.path.join(parent.data.path, item)
        fold = SniffFolder(path)
        node = treelib.Node(tag=fold.label, data=fold)
        tree.add_node(node, parent=parent.identifier)
        _add_children_nodes(tree, node)
    return


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


def _filter_path(path, node):
    abspath = os.path.abspath(os.path.expanduser(path))
    boo = False
    if node.data is not None:
        if node.data.path == abspath:
            boo = True
    return boo
