#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module defines the structure for filters objects (see FilterGeneral).

Various filters are then defined in submodules, by subclassing.

The subclass needs to define at least two things:
    * the attribute `_label` (string or unicode) which defines filter operation
    * the `func` method, that performs and outputs the boolean test

Some useful functions are defined here as well.
"""
from __future__ import print_function

import inspect  # introspection module
import numpy as np
import collections
try:
    from StringIO import StringIO  # python2
except ImportError:
    from io import StringIO  # python3


def bounded(arg, lower_bound=None, upper_bound=None):
    """Function that test whether argument is bounded.

    By convention, lower bound is included, upper bound is excluded. Thus it
    tests whether lower_bound <= arg < upper_bound. If arg is an iterable, the
    test must be satisfied for every element (hence the minimal value must be
    greater or equal to lower_bound and the maximal value must be lower than
    the upper bound).

    Parameters
    ----------
    arg : int or float
       quantity to be tested for bounds
    lower_bound : int or float (default None)
    upper_bound : int or float (default None)

    Returns
    -------
    Boolean
    """
    import collections
    # case : if we have an iterable, check bound on min, max values
    if isinstance(arg, collections.Iterable):
        if len(arg) == 0:  # empty sequence: False
            lower = False
            upper = False
        else:
            val_min = np.amin(arg)
            val_max = np.amax(arg)
            if lower_bound is not None:
                lower = val_min >= lower_bound
            else:
                lower = True
            if upper_bound is not None:
                upper = val_max < upper_bound
            else:
                upper = True
    # case : otherwise it should be a number
    else:
        if lower_bound is not None:
            lower = arg >= lower_bound
        else:
            lower = True
        if upper_bound is not None:
            upper = arg < upper_bound
        else:
            upper = True
    return (lower and upper)


def intersect(values, lower_bound=None, upper_bound=None):
    if not isinstance(values, list):
        values = [values, ]
    tmin = np.amin(values)
    tmax = np.amax(values)
    if lower_bound is not None:
        lower = tmin < upper_bound
    else:
        lower = True
    if upper_bound is not None:
        upper = tmax > lower_bound
    else:
        upper = True
    return lower and upper


def included(values, lower_bound=None, upper_bound=None):
    if not isinstance(values, list):
        values = [values, ]
    tmin = np.amin(values)
    tmax = np.amax(values)
    if lower_bound is not None:
        lower = tmin > lower_bound
    else:
        lower = True
    if upper_bound is not None:
        upper = tmax < upper_bound
    else:
        upper = True
    return lower and upper


class FilterError(Exception):
    "Superclass for errors while filtering"
    pass


class FilterTypeError(FilterError):
    "Error raised when one tries to combine different types of filters"
    pass


class FilterParamsError(FilterError):
    """Exception when parameters are not set for a given test."""
    pass


class FilterLabelError(FilterError):
    """Exception when label is not (properly) set."""
    pass


class FilterArgError(FilterError):
    """Exception raised when Argument format is not suitable."""
    pass


class FilterGeneral(object):
    """General class for filtering cell (i.e. tunacell.base.cell.Cell) instances.

    Important property is to make the instance callable, and define with a
    human readable label the action of filter.
    """
    _label = ''  # to be updated
    _type = None  # to be determined
    _obs = []  # declare observables that need to be computed prior to filter

    def __call__(self, *args):
        "Need `func` method -for each filter class- that performs boolean test"
        return self.func(*args)

    def _isattr(self, x):
        "Test whether x is NOt a method"
        return not inspect.ismethod(x)

    def __repr__(self):
        "Return string that can be called upon with built-in `eval` function"
        name = type(self).__name__
        chain = name + '('
        for name, val in inspect.getmembers(self, predicate=self._isattr):
            if name[0] != '_' and name != 'label':
                chain += '{}={}, '.format(name, repr(val))
        chain += ')'
        return chain

    def __str__(self):
        try:
            var = self._label
            if isinstance(var, str):
                return var
            else:
                raise FilterLabelError('No string representation')
        except AttributeError as ae:
            raise FilterLabelError(ae)

    def func(self, *args):
        """This is the boolean operation to define in specific Filter instances

        Default operation returns True.
        """
        return True

    @property
    def label(self):
        "Get label of applied filter(s)"
        if hasattr(self, '_label'):
            return self._label
        else:
            return None

    @label.setter
    def label(self, arg):
        "Set/update label"
        if not isinstance(arg, str):
            raise FilterArgError('The argument of the label setter \
                                 should be string/unicode')
        self._label = self._type + ', ' + arg
        return

    @property
    def obs(self):
        """Provides the list of hidden observables"""
        return self._obs


class FilterBoolean(FilterGeneral):
    """General class to implement Boolean operations between filters

    """

    _sequence = []  # seq. of filters with identical types (CELL, LINEAGE, ...)

    def __repr__(self):
        "specific repr since FilterGeneral.__repr__ does not work"
        name = type(self).__name__
        chain = name + '('
        for filt in self._sequence:
            chain += '{}, '.format(repr(filt))
        chain += ')'
        return chain

    def _update_obs(self):
        """Pulls up observable list"""
        ls = []
        for filt in self._sequence:
            for obs in filt._obs:
                if obs not in ls:
                    ls.append(obs)
        self._obs = ls
        return

    @property
    def label(self):
        "Get label of applied filter(s)"
        if hasattr(self, '_label'):
            return self._label
        else:
            return None

    @label.setter
    def label(self, arg):
        "Set/update label"
        if not isinstance(arg, str):
            raise FilterArgError('The argument of the label setter \
                                 should be string/unicode')
        self._label = arg
        return


class FilterTRUE(FilterBoolean):
    """Returns True for argument of any type, can be used in Boolean filters

    We need this artificial Filter to plug in defaults FilterSets.

    """

    def __init__(self):
        self._type = 'ANY'
        self.label = 'Always True'
        return

    def func(self, *args):
        return True


class FilterAND(FilterBoolean):
    """Defines boolean AND operation between same type filters.

    Parameters
    ----------
    filters : sequence of :class:`FilterGeneral` instances

    Returns
    -------
    :class:`FilterGeneral` instance
        will perform AND boolean operation between the various filters passed
        as arguments.

    Notes
    -----
    Defines a FilterTrue
    """
    def __init__(self, *filters):
        # register sequence of filters
        self._sequence = []
        types = []
        for filt in filters:
            if isinstance(filt, FilterAND):
                for ffilt in filt._sequence:
                    self._sequence.append(ffilt)
                    if ffilt._type not in types and ffilt._type != 'ANY':
                        types.append(ffilt._type)
            elif isinstance(filt, collections.Iterable):
                for ffilt in filt:
                    if not isinstance(ffilt, FilterGeneral):
                        msg = ('arg: {} is not a Filter'.format(ffilt))
                        raise FilterTypeError(msg)
                    self._sequence.append(ffilt)
                    if ffilt._type not in types and ffilt._type != 'ANY':
                        types.append(ffilt._type)
            elif isinstance(filt, FilterGeneral):
                self._sequence.append(filt)
                if filt._type not in types and filt._type != 'ANY':
                    types.append(filt._type)
            else:
                raise FilterTypeError('arg: {} is not a Filter'.format(filt))
        if len(types) > 1:
            raise FilterTypeError
        elif len(types) == 0:
            self._type = 'ANY'
        else:
            self._type = types[0]
        # set initial label
        label = ''
        # loop to complete label and check other things
        for index, filt in enumerate(self._sequence, start=0):
            if index == 0:
                indent = ' '*4
            else:
                indent = 'AND '
            for line in StringIO(filt.label):
                label += indent + line
            label += '\n'
        self.label = label.rstrip()
        # pulls up observable list from content
        self._update_obs()
        return

    def func(self, target):
        boo = True
        # when first False is found, no need to go further
        for filt in self._sequence:
            boo = filt(target)
            if not boo:
                break
        return boo


class FilterOR(FilterBoolean):
    """Defines boolean OR operation between same type filters.

    Parameters
    ----------
    filters : sequence of :class:`FilterGeneral` instances
    """

    def __init__(self, *filters):
        # register sequence of filters
        self._sequence = []
        types = []
        for filt in filters:
            if isinstance(filt, FilterOR):
                for ffilt in filt._sequence:
                    self._sequence.append(ffilt)
                    if ffilt._type not in types and ffilt._type != 'ANY':
                        types.append(ffilt._type)
            elif isinstance(filt, collections.Iterable):
                for ffilt in filt:
                    if not isinstance(ffilt, FilterGeneral):
                        msg = ('arg: {} is not a Filter'.format(ffilt))
                        raise FilterTypeError(msg)
                    self._sequence.append(ffilt)
                    if ffilt._type not in types and ffilt._type != 'ANY':
                        types.append(ffilt._type)
            elif isinstance(filt, FilterGeneral):
                self._sequence.append(filt)
                if filt._type not in types and filt._type != 'ANY':
                    types.append(filt._type)
            else:
                raise FilterTypeError('arg: {} is not a Filter'.format(filt))
        if len(types) > 1:
            raise FilterTypeError
        elif len(types) == 0:
            self._type = 'ANY'
        else:
            self._type = types[0]
        # set initial label
        label = ''
        # loop to complete label and check other things
        for index, filt in enumerate(self._sequence):
            if index == 0:
                indent = ' '*3
            else:
                indent = 'OR '
            for line in StringIO(filt.label):
                label += indent + line
            label += '\n'
        self.label = label.rstrip()
        # pulls up observable list from content
        self._update_obs()
        return

    def func(self, target):
        boo = False
        # when first True is found, no need to go further
        for filt in self._sequence:
            boo = filt(target)
            if boo:
                break
        return boo


class FilterNOT(FilterBoolean):
    """Defines boolean NOT operation on filter.

    Parameters
    ----------
    filter : :class:`FilterGeneral` instance
    """
    def __init__(self, filt):
        if not isinstance(filt, FilterGeneral):
            raise FilterTypeError('arg: {} is not a Filter'.format(filt))
        self._sequence = [filt, ]
        self._type = filt._type
        self.label = 'NOT [' + filt.label + ']'
        # pulls up observable list from content
        self._update_obs()
        return

    def func(self, target):
        filt, = self._sequence
        return not filt(target)


class FilterSet(object):
    """Collects filters of each type in a single object"""

    def __init__(self, label=None,
                 filtercell=FilterTRUE(),
                 filterlineage=FilterTRUE(),
                 filtertree=FilterTRUE(),
                 filtercontainer=FilterTRUE()):
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

        self._named_list = [('label', self.label),
                            ('filtercell', self.cell_filter),
                            ('filterlineage', self.lineage_filter),
                            ('filtertree', self.colony_filter),
                            ('filtercontainer', self.container_filter)]
        return

    def _collect_hidden_obs(self):
        """Returns the list of hidden observables (needed for computations)"""
        self._obs = []
        for filt in [self.cell_filter, self.lineage_filter, self.colony_filter, self.container_filter]:
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
        chain = name + '('
        for (name, item) in self._named_list:
            chain += '{}={}, '.format(name, repr(item))
        chain += ')'
        return chain

    def __str__(self):
        label = ''
        if self.label is not None:
            label += 'label: {}\n'.format(self.label)
        for filt in [self.cell_filter, self.lineage_filter, self.colony_filter,
                     self.container_filter]:
            if filt._type == 'ANY':
                continue
            for index, line in enumerate(StringIO(str(filt))):
                if index == 0:
                    label += '* '
                else:
                    label += '  '
                label += line
            label += '\n'
        return label.rstrip()
    