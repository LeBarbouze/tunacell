#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module defines core classes storing/reading statistical analysis of the
dynamics of a single observable.

UnivariateConditioned(obs, applied_filter=None)
    Object that stores data for a given condition (a condition is defined by
    the 'applied_filter' Filter instance given as parameter).

Univariate(obs, cset=[], parser=None, indexify=None)
    Object that stores data for a given set of conditions as a dictionary of
    UnivariateConditioned instances, whose keys are each filter string
    representation.

StationaryUnivariateConditioned(obs, applied_filter=None)
    Object that stores data under stationary hypothesis for a given condition
    (a condition is defined by the 'applied_filter' Filter instance given as
    parameter).

StationaryUnivariate(obs, cset=[], parser=None)
    Object that stores data under stationary hypothesis for a given set of
    conditions, as a dictionary of StationaryUnivariateConditioned instances,
    whose keys are each filter string representation.
"""
from __future__ import print_function

import os

import numpy as np
import pandas as pd
import warnings
from tuna.io import text

from tuna.stats.utils import Regions, CompuParams


class UnivariateConditioned(object):
    """Stores statistics of dynamics of a single observable, single condition.

    This is where one-point and two-point functions are stored. Such objects
    usually live as entries in a global Univariate instance.

    Parameters
    ----------
    obs : Observable instance
    applied_filter : Filter instance (that defines a condition)
    indexify : :class:`Indexify` instance
        used to bin time values

    Attributes
    ----------
    time : 1d array of floats
        (ordered) sequence of times at which data is collected (binned)
    count_one : 1d array of ints
        number of samples for one point function, associated to corresponding
        time array.
    average : 1d array of floats
        sample average associated to time array
    var : 1d array of floats
        sample variance associated to time array
    std : 1d array of floats
        sample standard deviation (square root of the variance)
    count_two : 2d array of ints
        matrix of two-point counts (number of trajectory connecting times i and
        j, where i and j are row, column indices to corresponding time in time
        array)
    autocorr : 2d array of floats
        matrix of sample autocovariance between times i and j, associated to
        corresponding entries in time array.

    See also
    --------
    * definition of tuna.stats.compute.set_dynamics()
    """

    def __init__(self, univariate, applied_filter=None):
        self.univariate = univariate  # parent object
        self.applied_filter = applied_filter
        if applied_filter is not None:
            self.condition = repr(applied_filter)  # can be called later...
        else:
            self.condition = 'master'
        # define attributes that store data
        self.onepoint = None  # Numpy structured array (time, count, av, sd)
        self.count_two = None  # 2d array
        self.autocorr = None  # 2d array

        self._keys = ('array', 'count_two', 'autocorr')
        self.stationary = None  # StationaryUnivariateConditioned instance
        self.two = {}  # stores cross-correlation results: key=str(other obs)
        self.two_stationary = {}
        return

    def bind(self, array, count_two, autocorr):
        self.onepoint = array  # numpy array: time, count, average, std_dev
        self.count_two = count_two
        self.autocorr = autocorr
        return

    @property
    def time(self):
        if self.onepoint is not None:
            return self.onepoint['time']
        else:
            return None

    @property
    def count_one(self):
        if self.onepoint is not None:
            return self.onepoint['count']
        else:
            return None

    @property
    def average(self):
        if self.onepoint is not None:
            return self.onepoint['average']
        else:
            return None

    @property
    def var(self):
        if self.autocorr is not None:
            return np.diagonal(self.autocorr)
        else:
            return None

    @property
    def std(self):
        if self.autocorr is not None:
            return np.sqrt(self.var)
        else:
            return None

    def _get_obs_path(self, user_root=None, write=False):
        obs = self.univariate.obs
        exp = self.univariate.parser.experiment
        fset = self.univariate.parser.fset
        analysis_path = text.get_analysis_path(exp, user_abspath=user_root,
                                               write=write)
        res = text.get_filter_path(analysis_path, fset, write=write)
        index_filter, filter_path = res
        res = text.get_condition_path(filter_path, self.applied_filter,
                                      write=write)
        index_condition, condition_path = res
        obs_path = text.get_observable_path(condition_path, obs, write=write)
        return obs_path

    def write_text(self, path=None):
        """Write arrays to files

        Parameters
        ----------
        path : str (default None)
            analysis folder path under which filterset->condition->obs
            leave to None to canonical analysis path under the experiment
            analysis folder
        """
        # check path and write files
        obs_path = self._get_obs_path(user_root=path, write=True)
        ffmt = '%.8e'  # floating point numbers
        ifmt = '%d'  # integers
        item_path = os.path.join(obs_path, 'onepoint.tsv')
        names = self.onepoint.dtype.names
        header = '\t'.join(names)
        fmt = [ifmt if 'count' in n_ else ffmt for n_ in names]
        np.savetxt(item_path, self.onepoint, fmt=fmt,
                   delimiter='\t', comments='', header=header)

        for key in ['count_two', 'autocorr']:
            array = self[key]
            item_path = os.path.join(obs_path, key + '.tsv')
            if 'count' in key:
                fmt = ifmt
            else:
                fmt = ffmt
            np.savetxt(item_path, array, fmt=fmt, delimiter='\t')
        return

    def read_text(self, path=None):
        """Initialize object by reading text output."""
        # check path and write files
        obs_path = self._get_obs_path(user_root=path, write=False)
        # read
        item_path = os.path.join(obs_path, 'onepoint.tsv')
        if not os.path.exists(item_path):
            raise text.MissingFileError(item_path)
        array = np.genfromtxt(item_path, delimiter='\t',
                              dtype=(float, int, float, float), names=True)
        self.onepoint = array
        for key in ['count_two', 'autocorr']:
            if 'count' in key:
                dtype = int
            else:
                dtype = float
            item_path = os.path.join(obs_path, key + '.tsv')
            if not os.path.exists(item_path):
                raise text.MissingFileError(item_path)
            array = np.genfromtxt(item_path, delimiter='\t', dtype=dtype)
            self[key] = array
        return

    def __setitem__(self, key, val):
        if key not in self._keys:
            msg = 'key "{}" is not valid.'.format(key)
            msg += '\nChoose a key in following list:'
            for key in self._keys:
                msg += '\n- {}'.format(key)
            warnings.warn(msg)
            val = None
        else:
            self.__setattr__(key, val)
        return

    def __getitem__(self, key):
        if key not in self._keys:
            msg = 'key "{}" is not valid.'.format(key)
            msg += '\nChoose a key in following list:'
            for key in self._keys:
                msg += '\n- {}'.format(key)
            warnings.warn(msg)
            val = None
        else:
            val = self.__getattribute__(key)
        return val

    def info(self):
        """Prints info on stored data."""
        if self.time is None:
            msg = (' NO VALUES')
            return msg
        msg = (' One-point arrays' + '\n'
               '  Times:' + '\n'
               '   {}'.format(self.time) + '\n'
               '   ({} values)'.format(len(self.time)) + '\n'
               '   min value {}'.format(np.amin(self.time)) + '\n'
               '   max value {}'.format(np.amax(self.time)) + '\n'
               '  Sample counts:' + '\n'
               '   {}'.format(self.count_one) + '\n'
               '   ({} values)'.format(len(self.count_one)) + '\n'
               '   min value {}'.format(np.amin(self.count_one)) + '\n'
               '   max value {}'.format(np.amax(self.count_one)) + '\n'
               '   average   {}'.format(np.mean(self.count_one)) + '\n'
               '  Sample average:' + '\n'
               '   {}'.format(self.average) + '\n'
               '   ({} values)'.format(len(self.average)) + '\n'
               '   min value {}'.format(np.amin(self.average)) + '\n'
               '   max value {}'.format(np.amax(self.average)) + '\n'
               '   average   {}'.format(np.mean(self.average)) + '\n'
               ' Two-point matrices:' + '\n'
               '   Sample counts:' + '\n'
               '    shape    {}'.format(self.count_two.shape) + '\n'
               )
        return msg
#
#    def _compute_naive_stationary(self, tmin=None, tmax=None):
#        """Computes stationary auto-correlation between tmin and tmax.
#
#        Beware that this computation is of poor accuracy, as it relies solely
#        on the dynamical autocorrelation matrix, the values of which are
#        proned to sampling errors. A better computational tool is provided
#        in :mod:`stats.compute`, as :func:`set_stationary_autocorrelation`,
#        or more accessibly through:mod:`stats.api` with
#        :func:`under compute_observable_stationary`.
#        """
#        dts, cts, res = get_stat_from_dynamics(self, tmin, tmax)
#        sindexify = deepcopy(self.indexify)
#        sindexify.offset = 0.  # to offseting for time differences
#
#        array = np.zeros(len(dts), dtype=[('time_interval', 'f8'),
#                                          ('count', 'u4'),
#                                          ('auto-correlation', 'f8')])
#        array['time_interval'] = dts
#        array['count'] = cts
#        array['auto-correlation'] = res
#
#        SOC = StationaryUnivariateConditioned  # short alias
#
#        stat = SOC(self.obs, applied_filter=self._applied_filter,
#                   indexify=sindexify, tmin=tmin, tmax=tmax,
#                   array=array)
#        self.naive_stationary = stat
#
#        return stat


class UnivariateIOError(IOError):
    pass


class Univariate(object):
    """Object that stores the dynamical analysis of a given observable.

    This object is created when the function compute_observable_dynamics() is
    called. It stores single condition UnivariateConditioned instances,
    for all conditions provided in cset argument, in addition to 'master' (no
    condition).

    Parameters
    ----------
    obs : :class:`Observable` instance
    cset : list of :class:`FilterSet` instances (default [])
        default: only 'master' is computed
    parser : :class:`Parser` instance
        Defines which experiment, and how to parse it (with its filter set
        fset attribute)
    region : object
        with tmin, tmax, name attributes; first two are used to set boundaries
        for absolute time data (accept only acquisitions within those bounds),
        'name' is used for storage
    period : float

    Attributes
    ----------
    _items : dictionary
        keys are condition labels, values are
        class:`UnivariateConditioned` instances.
    """

    def __init__(self, obs, cset=[], parser=None, region=None, eval_times=None):
        self.obs = obs
        self.parser = parser
        self.exp = parser.experiment
        self.cset = cset
        self.region = region
        self.eval_times = eval_times
        # create as many nodes as there are conditions in cset
        self._items = {}
        self._condition_labels = []
        # alias
        Unic = UnivariateConditioned
        # Instantiate a SingleObsConditioned, without applied_filter: master
        master = Unic(self, applied_filter=None)
        # store it as 'master'
        self._items['master'] = master
        self._condition_labels.append('master')
        for cdt in cset:
            cdt_repr = '{}'.format(repr(cdt))
            self._condition_labels.append(cdt_repr)
            self._items[cdt_repr] = Unic(self, cdt)
        return

    def __getitem__(self, key):
        return self._items[key]

    def __setitem__(self, key, val):
        self._items[key] = val
        return

    @property
    def master(self):
        """There's always a master (no condition)"""
        return self['master']

    def __str__(self):
        msg = ('--------------------------' + '\n'
               'STATISTICS OF THE DYNAMICS' + '\n'
               '--------------------------' + '\n'
               '' + '\n'
               'OBSERVABLE: {}'.format(self.obs.label) + '\n'
               'CONDITIONS:' + '\n'
               '  * `master` (no condition)' + '\n'
               )
        for condition_label in self._condition_labels[1:]:
            msg += '  * `{}`'.format(condition_label) + '\n'
        msg += 'BINNING:\n{}'.format(self.indexify)
        msg += '\n'
        msg += 'Showing info for `master`:' + '\n'
        msg += str(self._items['master'].info())
        return msg

#    def compute_naive_stationary(self, tmin=None, tmax=None):
#        """Compute stationary auto-correlation for all conditions.
#        """
#        so = StationaryUnivariate(self.obs, cset=self.cset, parser=self.parser,
#                                  tmin=tmin, tmax=tmax)
#        for key in self._condition_labels:
#            stat = self[key]._compute_naive_stationary(tmin, tmax)
#            so[key] = stat
#        self.naive_stationary = so
#        return

    def export_text(self, analysis_folder=None):
        """Export results to text files.

        Parameters
        ----------
        analysis_folder : str (default None)
            Path to the analysis folder; default is 'analysis' subfolder in
            experiment folder
        """
#        analysis_path = text.get_analysis_path(self.parser.experiment,
#                                               user_abspath=analysis_folder,
#                                               write=True)
#        # write Indexify associated to obs
#        obs = self.obs
#        fset = self.parser.fset
#        index_filter, filter_path = text.get_filter_path(analysis_path, fset,
#                                                         write=True)
#        basename = 'indexify_' + obs.label + '.txt'
#        text_file = os.path.join(filter_path, basename)
#        with open(text_file, 'w') as f:
#            f.write(repr(self.indexify))
        # write each condition
        for key, val in self._items.items():
            val.write_text(analysis_folder)
        return

    def import_from_text(self, analysis_folder=None):
        """Set instance from text files.

        This reader needs that Observable, Parser, Conditions are defined
        and must match text folders.

        Raises
        ------
        UnivariateIOError
            when any of the folder is not found, or indexify does not match

        """
#        analysis_path = text.get_analysis_path(self.parser.experiment,
#                                               user_abspath=analysis_folder,
#                                               write=False)
#        if not os.path.exists(analysis_path):
#            raise UnivariateIOError('No analysis folder')
#        obs = self.obs
#        fset = self.parser.fset
#        try:
#            _, filter_path = text.get_filter_path(analysis_path, fset,
#                                                  write=False)
#        except text.MissingFolderError as missing:
#            raise UnivariateIOError(missing)
#        if not os.path.exists(filter_path):
#            raise UnivariateIOError('Missing folder {}'.format(filter_path))
#        # reading Indexify file
#        basename = 'indexify_' + obs.label + '.txt'
#        text_file = os.path.join(filter_path, basename)
#        try:
#            with open(text_file, 'r') as f:
#                rep = f.readline()
#        except IOError:  # no file, repr of indexify is empty
#            rep = ''
#        except text.MissingFolderError as missing:
#            raise UnivariateIOError(missing)
#        if rep.rstrip() != repr(self.indexify):
#            raise UnivariateIOError('Problem at Indexify')
        # read each condition
        try:
            for key, val in self._items.items():
                val.read_text(analysis_folder)
        except (text.MissingFileError, text.MissingFolderError) as missing:
            raise UnivariateIOError(missing)
        return


# %% STATIONARY AUTOCORRELATION
class StationaryUnivariateConditioned(object):
    """Stores stationary analysis for a given observable, with a condition.

    Similar structure than ObservableConditioned, although all values are
    recorded in a single Numpy array, with 3 columns.

    Parameters
    ----------
    obs : :class:`Observable` instance
    label : str (default None)
    applied_filter : :class:`FilterSet` instance
    indexify : :class:`Indexify` instance
    tmin : float (default None)
    tmax : float (default None)
    adjust_mean : str {'global', 'local'}
    array : 3-columns structured array (default None)
        where results will be stored, set to None by default.
        Column names are ['time_interval', 'count', 'auto_correlation']

    Attributes
    ----------
    self.array : 3 columns Numpy array (time_interval, counts, autocorr)

    See also
    --------
    StationaryUnivariate:
        a collection of StationaryUnivariateConditioned
        instances where export/import are defined.
    """
    def __init__(self, statunivariate, applied_filter=None, array=None):
        self.statunivariate = statunivariate
        self.basename = 'stationary'
        # add region label
        self.basename += '_' + self.statunivariate.region.name
        # add computation options
        self.basename += '_' + self.statunivariate.options.as_string_code()
        # condition
        self.applied_filter = applied_filter
        if applied_filter is not None:
            self.condition = repr(applied_filter)
        else:
            self.condition = 'master'
        self.array = array  # should be a 3 column array of dtype:
        # [('time_interval', 'f8'),('count', 'i4'), ('autocorrelation', 'f8')]
        return

    @property
    def time(self):
        if self.array is not None:
            return self.array['time_interval']
        else:
            return None

    @property
    def count(self):
        if self.array is not None:
            return self.array['count']
        else:
            return None

    @property
    def autocorr(self):
        if self.array is not None:
            return self.array['auto_correlation']
        else:
            return None

    def write_text(self, path='.'):
        """Write array to file."""
        # get observable path that should exist already
        cdt_univ = self.statunivariate.univariate[self.condition]
        obs_path = cdt_univ._get_obs_path(user_root=path, write=False)
        # otherwise a text.MissingFolder will be raised here
        if self.array is None:
            print('Nothing to write')
            return
        ffmt = '%.8e'  # floating point numbers
        ifmt = '%d'  # integers
        item_path = os.path.join(obs_path, self.basename + '.tsv')
        names = self.array.dtype.names
        header = '\t'.join(names)
        fmt = [ifmt if 'count' in n_ else ffmt for n_ in names]
        np.savetxt(item_path, self.array, fmt=fmt,
                   delimiter='\t', comments='', header=header)
        return

    def read_text(self, path='.'):
        """Initialize object by reading text output."""
        # get observable path that should exist already
        cdt_univ = self.statunivariate.univariate[self.condition]
        obs_path = cdt_univ._get_obs_path(user_root=path, write=False)
        item_path = os.path.join(obs_path, self.basename + '.tsv')
        if not os.path.exists(item_path):
            raise text.MissingFileError(item_path)
        arr = np.genfromtxt(item_path, delimiter='\t',
                            dtype=(float, int, float, float),
                            names=True)
        self.array = arr
        return

    def info(self):
        """Return string giving info on stored data."""
        if self.array is None:
            msg = (' NO VALUES')
            return msg
        times = self.array['time_interval']
        counts = self.array['count']
        autocorrs = self.array['auto_correlation']
        msg = ('  Time intervals:' + '\n'
               '   {}'.format(times) + '\n'
               '   ({} values)'.format(len(times)) + '\n'
               '   min value {}'.format(np.amin(times)) + '\n'
               '   max value {}'.format(np.amax(times)) + '\n'
               '  Sample counts:' + '\n'
               '   {}'.format(counts) + '\n'
               '   min value {}'.format(np.amin(counts)) + '\n'
               '   max value {}'.format(np.amax(counts)) + '\n'
               '   average   {}'.format(np.mean(counts)) + '\n'
               '  Autocorrelation values:' + '\n'
               '   {}'.format(autocorrs) + '\n'
               '   ({} values)'.format(len(autocorrs)) + '\n'
               '   min value {}'.format(np.amin(autocorrs)) + '\n'
               '   max value {}'.format(np.amax(autocorrs)) + '\n'
               '   average   {}'.format(np.mean(autocorrs)) + '\n'
               )
        return msg


class StationaryUnivariateIOError(IOError):
    pass


class StationaryUnivariate(object):
    """Object that stores the stationary analysis of a given observable.

    Similar structure than for Univariate.

    This object is created when the function compute_observable_stationary() is
    called. It stores StationaryUnivariateConditioned instances,
    for all conditions provided in cset argument, in addition to 'master' (no
    condition).

    Parameters
    ----------
    univariate : :class:`Univariate` instance
        from which one needs to compute stationary autocorrelation function
    region : pandas.Series
        name is the label of the region, and two items under index [tmin, tmax]
    options : :class:`tuna.stats.utils.CompuParams` instance
        define how to substract mean value and how to accept segments

    Attributes
    ----------
    univariate : :class:`Univariate` instance
        from which one needs to compute stationary autocorrelation function
    region : pandas.Series
        name is the label of the region, and two items under index [tmin, tmax]
    tmin : float (default None)
        minimum time (included) allowed in collecting timeseries
    tmax : float (default None)
        maximum time (excluded) allowed in collecting timeseries
    options : :class:`tuna.stats.utils.CompuParams` instance
        define how to substract mean value and how to accept segments
    adjust_mean : str {'global', 'local'}
        whether to use local average estimate (average estimated at single
        time points), or global average estimate (pool of all estimated values
        between tmin and tmax)
    disjoint : bool {True, False}
    _items : dictionary
        keys are condition labels, values are StationaryUnivariateConditioned
        instances.

    See also
    --------
    compute_naive_stationary_observable
        function to compute autocorrelation function with stationary hypothesis
    """

    def __init__(self, univariate, region=None, options=None):
        self.univariate = univariate
        if region is not None:
            self.region = region
        else:
            regs = Regions(univariate.exp)
            self.region = regs.get('ALL')  # get default region
        if options is not None:
            self.options = options
        else:
            self.options = CompuParams()
        self.obs = univariate.obs
        self.parser = univariate.parser
        self.cset = univariate.cset
        self.label = self.region.name
        self.tmin = self.region.tmin
        self.tmax = self.region.tmax
        self.adjust_mean = self.options.adjust_mean
        self.disjoint = self.options.disjoint
        # pandas.DataFrame to store values in table
        self.dataframe = None
        # create as many nodes as there are conditions in cset
        self._items = {}
        self._condition_labels = []
        # alias
        SUnic = StationaryUnivariateConditioned
        # add item for no conditions, 'master'
        self._items['master'] = SUnic(self, applied_filter=None, array=None)
        self._condition_labels.append('master')
        for cdt in univariate.cset:
            cdt_repr = '{}'.format(repr(cdt))
            self._condition_labels.append(cdt_repr)
            self._items[cdt_repr] = SUnic(self, applied_filter=cdt, array=None)
        return

    def __getitem__(self, key):
        return self._items[key]

    def __setitem__(self, key, val):
        self._items[key] = val
        return

    @property
    def master(self):
        """There's always a master (no condition)"""
        return self['master']

    def __str__(self):
        msg = ('-----------------------------------------------' + '\n'
               'STATISTICS OF THE DYNAMICS: STATIONARY ANALYSIS' + '\n'
               '-----------------------------------------------' + '\n'
               '' + '\n'
               'OBSERVABLE: {}'.format(self.obs.label) + '\n'
               'CONDITIONS:' + '\n'
               '  * `master` (no condition)' + '\n'
               )
        for condition_label in self._condition_labels[1:]:
            msg += '  * `{}`'.format(condition_label) + '\n'
        msg += '\n'
        msg += 'Showing info for `master`:' + '\n'
        msg += str(self._items['master'])
        return msg

    def export_text(self, analysis_folder=None):
        # write each condition
        try:
            for key, val in self._items.items():
                val.write_text(analysis_folder)
        # when not possible it means single object has not been exported yet
        except text.MissingFolderError:
            self.single.export_text(analysis_folder)
            for key, val in self._items.items():
                val.write_text(analysis_folder)
        # export dataframe as csv file
        if self.dataframe is not None:
            exp = self.univariate.parser.experiment
            fset = self.univariate.parser.fset
            analysis_path = text.get_analysis_path(exp, user_abspath=analysis_folder,
                                                   write=True)
            res = text.get_filter_path(analysis_path, fset, write=True)
            index_filter, filter_path = res
            basename = 'data_{}_{}'.format(self.label, self.obs.name)
            text_file = os.path.join(filter_path, basename + '.csv')
            self.dataframe.to_csv(text_file, index=False)
        return

    def import_from_text(self, analysis_folder=None):
        try:
            for key, val in self._items.items():
                val.read_text(analysis_folder)
            exp = self.univariate.parser.experiment
            fset = self.univariate.parser.fset
            analysis_path = text.get_analysis_path(exp, user_abspath=analysis_folder,
                                                   write=False)
            res = text.get_filter_path(analysis_path, fset, write=False)
            index_filter, filter_path = res
            basename = 'data_{}_{}'.format(self.label, self.obs.name)
            text_file = os.path.join(filter_path, basename + '.csv')
            self.dataframe = pd.read_csv(text_file, index_col=False)
        except (text.MissingFileError, text.MissingFolderError):
            raise StationaryUnivariateIOError
        return
