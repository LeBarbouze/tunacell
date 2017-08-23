#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module provides API class definition to define observables.

Classes
-------

* :class:`Observable`: main object to define observables.

"""
from __future__ import print_function

import warnings
import re
import inspect


# DEPRECATED
_re_codestring = 'T([a-z])(\d*[\.,]*\d*)M([a-z\-]+)J(\d+)'


def _is_valid_codestring(codestring):
    chain = re.compile(_re_codestring)
    m = chain.match(codestring)
    if m:
        return True
    else:
        return False


class ObservableError(Exception):
    pass


class ObservableStringError(ObservableError):
    pass


class Observable(object):
    """Defines how to retrieve observables.

    Parameters
    ----------
    label : str
        user name for this observable (can be one of the raw observable)
    raw : str (default None)
        raw name of the observable: must be a column name of raw data, i.e.
        first element of one entry of Experiment.datatype
    differentiate : boolean (default False)
        whether to differentiate raw observable
    scale : str {'linear', 'log'}
        expected scaling form as a function of time
        (used for extrapolating values at boundaries, including in the
        local_fit procedure)
    local_fit : boolean (default False)
        whether to perform local fit procedure
    time_window : float (default 20.)
        Time window over which local fit procedure is applied
        (used only when local_fit is activated)
    join_points : int (default 3)
        number of points over which extrapolation procedure is led by
        fitting linearly [the scale of] the observable w.r.t time
    mode : str {'dynamics', 'birth', 'division', 'net_increase', 'rate',
        'average'}
        mode used to retrieve data:
            * 'dynamics': all timepoints are retrieved
            * 'birth': only birth value is retrieved
            * 'division': only division value is retrieved
            * 'net-increase-additive': difference between division value
               and birth value
            * 'net-increase-multiplicative': ratio between division value
               and birth value
            * 'rate': rate of linear fit of [scale of] observable
            * 'average': average of observable over cell cycle
    timing : str {'t', 'b', 'd', 'm', 'g'}
        set the time at which cell cycle observable is associated:
            * 't' : time-lapse timing (associated to mode 'dynamics')
            * 'b' : cell cycle birth time
            * 'd' : cell cycle division time
            * 'm' : cell cycle midpoint (half time)
            * 'g' : cell cycle generation index
    tref : float
        when timing is set to 'g', sets the 0th generation to the cell
        that bounds this reference time
    """

    def __init__(self, label,
                 raw=None, differentiate=False, scale='linear',
                 local_fit=False, time_window=0., join_points=3,
                 mode='dynamics', timing='t', tref=None):
        self._attr_names = ['name',
                            'raw',
                            'scale',
                            'differentiate',
                            'local_fit',
                            'time_window',
                            'join_points',
                            'mode',
                            'timing',
                            'tref']

        self._label = label
        if raw is None:
            self.raw = label
        else:
            self.raw = raw
        self.differentiate = differentiate
        self.scale = scale
        self.local_fit = local_fit
        self.time_window = time_window
        self.join_points = join_points
        self.mode = mode
        self.timing = timing
        self.tref = tref
        return

    @property
    def label(self):
        return self._label

    # DEPRECATED
    def old_label(self):
        """Label is outputing a string representation suitable for file label.

        This method creates a string label that specifies each parameter to
        re-construct the Observable. The output string is
        (kind-of) human readable. More importantly, it is suitable
        to be a filename (only alphanumeric caracters, and underscores), and in
        fact serves to name directories in analysis folder.

        Note
        ----
        :func:`__repr__` : returns another string representation, that can be
        called by the built-in :func:`eval()`, to instantiate a new object
        with identical functional parameters.
        """
        msg = ''
        # timing is in between T flags
        if self.tref is not None:
            stref = '{}'.format(self.tref)
        else:
            stref = ''
        msg += 'T' + self.timing + stref
        # mode is in between M flags
        msg += 'M' + self.mode
        msg += 'J' + '{}'.format(self.join_points)
        msg += '_'
        if self.differentiate:
            msg += 'dot_'
        if self.local_fit:
            msg += 'W{:.2f}_'.format(self.time_window)
        if self.scale == 'log':
            msg += 'log_'
        msg += self.raw
        return msg

    # DEPRECATED
    def load_from_string(self, codestring):
        """Set Observable instance from string code created by `label` method

        Parameters
        ----------
        codestring : str
            must follow some rules for parsing
        """
        # set options to default and update if found
        self.mode = 'dynamics'
        self.timing = 't'
        self.local_fit = False
        self.time_window = 0.
        self.join_points = 3  # default
        self.differentiate = False
        self.scale = 'linear'
        items = codestring.split('_')
        self.raw = items[-1]  # last item is always raw observable label
        if self.raw == 'none':
            msg = ("'raw' is set to 'none'.\n"
                   "Update to a valid column name of your experiment.")
            warnings.warn(msg)

        # test whether codestring is valid: must have T and M flags
        chain = re.compile(_re_codestring)
        m = chain.match(codestring)
        if m:
            timing, stref, mode, sjoin = m.groups()
            self.timing = timing
            if stref:
                self.tref = float(stref.replace(',', '.'))  # if decimal is ,
            else:
                self.tref = None
            self.mode = mode
            self.join_points = int(sjoin)
        else:
            raise ObservableStringError('Not a valid codestring')

        # try to check whether local fit is performed and its parameters
        pfit = re.compile('W(\d*[.,]*\d*)')

        for item in items[:-1]:
            # local_fit?
            m = pfit.search(item)
            if m is not None:
                stime_window, = m.groups()
                # check that tw_str is not empty
                if stime_window:
                    self.time_window = float(stime_window.replace(',', '.'))
                    self.local_fit = True
            # log scale
            if item == 'log':
                self.scale = 'log'
            if item == 'dot':
                self.differentiate = True
        return

    def as_string_table(self):
        """Human readable output as a table.
        """
        msg = ''
        column_size = 14

        def formatting(items):
                out = ' | '.join(['{}'.format(items[0]).rjust(column_size),
                                  '{}'.format(items[1]).ljust(column_size)])
                return out

        msg += formatting(['parameter', 'value'])
        msg += '\n' + formatting(['----', '----'])
        for attr in self._attr_names:
            val = self.__getattribute__(attr)
            msg += '\n' + formatting([attr, val])
        msg += '\n'
        return msg

    def as_latex_string(self):
        """Export as LaTeX string.
        """
        output = r'$'
        if self.differentiate:
            output += '\\frac{\\mathrm{d}}{\\mathrm{d}t}'
            if self.scale == 'log':
                output += '\\log\\left('  # parenthesis started
        variable_name = '{}'.format(self.raw)
        output += '\\mathrm{{ {} }}'.format(variable_name.replace('_', '\, '))
        if self.timing == 't':
            output += '(t)'
        elif self.timing == 'b':
            output += '\\left( t_{\\mathrm{birth}} \\right)'
        elif self.timing == 'd':
            output += '\\left( t_{\\mathrm{div}} \\right)'
        elif self.timing == 'm':
            output += ('\\left( \\frac{t_{\\mathrm{birth}} + t_{\\mathrm{div}}}'
                       '{2} \\right)')
        elif self.timing == 'g':
            output += '\\left( \\mathrm{generation} \\right)'
        if self.differentiate and self.scale == 'log':
            output += '\\right)'  # parenthesis closed
        if self.mode != 'dynamics':
            output += '_{{\mathrm{{ {} }} }}'.format(self.mode)
        if self.local_fit:
            output += '\\ [window: {}]'.format(self.time_window)
        output += '$'
        return output

    def __str__(self):
        return self.label

    def __repr__(self):
        name = type(self).__name__
        chain = name + '('
        for attr in self._attr_names:
            val = self.__getattribute__(attr)
            chain += '{}={}, '.format(attr, repr(val))
        chain += ')'
        return chain


class FunctionalObservable(object):
    """Combination of :class:`Observable` instances"""

    def __init__(self, name, f, observables):
        self.name = name
        if not callable(f):
            raise ValueError('f must be callable')
        self.f = f
        argspec = inspect.getargspec(f)
        self.observables = observables
        if len(observables) != len(argspec.args):
            msg = ('length of observable list must match number of arguments of f ')
            raise ValueError(msg)
        return


if __name__ == '__main__':
    length = Observable('length', raw='length')
    width = Observable('width', raw='width')
    def volume(x, y):
        return x * y**2
    combo = FunctionalObservable('volume', volume, [length, width])
    