#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module defines core classes storing/reading statistical analysis of the
dynamics of a single observable.

UnivariateConditioned(obs, applied_filter=None)
    Object that stores data for a given condition (a condition is defined by
    the 'applied_filter' Filter instance given as parameter).

Univariate(obs, cset=[], exp=None, indexify=None)
    Object that stores data for a given set of conditions as a dictionary of
    UnivariateConditioned instances, whose keys are each filter string
    representation.

StationaryUnivariateConditioned(obs, applied_filter=None)
    Object that stores data under stationary hypothesis for a given condition
    (a condition is defined by the 'applied_filter' Filter instance given as
    parameter).

StationaryUnivariate(obs, cset=[], exp=None)
    Object that stores data under stationary hypothesis for a given set of
    conditions, as a dictionary of StationaryUnivariateConditioned instances,
    whose keys are each filter string representation.
"""
from __future__ import print_function

import os
import logging

import numpy as np
import pandas as pd
import warnings
from tunacell.io import analysis
from tunacell.base.observable import ObservableNameError
from tunacell.stats.utils import Region, Regions, CompuParams, _dtype_converter


logger = logging.getLogger(__name__)


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
    * definition of tunacell.stats.compute.set_dynamics()
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

    def _get_path(self, user_root=None, write=False):
        """Get condition path"""
        obs_path = self.univariate._get_obs_path(user_root=user_root, write=write)
        res = analysis.get_condition_path(obs_path, self.applied_filter, write=write)
        index_condition, condition_path = res
        return condition_path

    def bind(self, array, count_two, autocorr):
        self.onepoint = array  # numpy array: time, counts, average, std_dev
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
            return self.onepoint['counts']
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

    def write_text(self, path=None):
        """Write arrays to files

        Parameters
        ----------
        path : str (default None)
            analysis folder path under which filterset->condition->obs
            leave to None to canonical analysis path under the experiment
            analysis folder
        """
        ext = '.tsv'
        add_name = '_' + self.univariate.region.name
        condition_path = self._get_path(user_root=path, write=True)
        ffmt = '%.8e'  # floating point numbers
        ifmt = '%d'  # integers
        bname = 'onepoint' + add_name + ext
        fname = os.path.join(condition_path, bname)
        names = self.onepoint.dtype.names
        header = '\t'.join(names)
        fmt = [ifmt if 'count' in n_ else ffmt for n_ in names]
        np.savetxt(fname, self.onepoint, fmt=fmt,
                   delimiter='\t', comments='', header=header)

        for key in ['count_two', 'autocorr']:
            array = self[key]
            bname = key + add_name + ext
            fname = os.path.join(condition_path, bname)
            if 'count' in key:
                fmt = ifmt
            else:
                fmt = ffmt
            np.savetxt(fname, array, fmt=fmt, delimiter='\t')
        return

    def read_text(self, path=None):
        """Initialize object by reading text output."""
        ext = '.tsv'
        add_name = '_' + self.univariate.region.name
        condition_path = self._get_path(user_root=path, write=False)
        # read
        bname = 'onepoint' + add_name + ext
        fname = os.path.join(condition_path, bname)
        if not os.path.exists(fname):
            raise analysis.MissingFileError(fname)
        array = np.genfromtxt(fname, delimiter='\t',
                              dtype=(float, int, float, float), names=True)
        self.onepoint = array
        for key in ['count_two', 'autocorr']:
            if 'count' in key:
                dtype = int
            else:
                dtype = float
            bname = key + add_name + ext
            fname = os.path.join(condition_path, bname)
            if not os.path.exists(fname):
                raise analysis.MissingFileError(fname)
            array = np.genfromtxt(fname, delimiter='\t', dtype=dtype)
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

    def _onepoint_as_dataframe(self):
        return pd.DataFrame(self.onepoint)

    def _twopoint_as_dataframe(self):
        # unroll
        n_items = len(self.time)
        row_times = np.concatenate([np.array(n_items * [self.time[i], ]) for i in range(n_items)])
        col_times = np.concatenate([self.time for item in range(n_items)])
        counts = self.count_two.flatten()
        autocov = self.autocorr.flatten()
        names = ['time-row', 'time-col', 'counts', 'autocovariance']
        data = {'time-row': row_times,
                'time-col': col_times,
                'counts': counts,
                'autocovariance': autocov,
                }
        return pd.DataFrame(data, columns=names)

    def display_onepoint(self, item_max=10):
        print(self._onepoint_as_dataframe()[:item_max])
        return

    def display_twopoint(self, item_max=10, sampling=False):
        df = self._twopoint_as_dataframe()
        if sampling:
            show = df.sample(item_max).sort_index()
        else:
            show = df[:item_max]
        print(show)
        return


class UnivariateIOError(IOError):
    pass


class UnivariateInitError(ValueError):
    pass


class Univariate(object):
    """Object that stores the dynamical analysis of a given observable.

    This object is created when the function compute_observable_dynamics() is
    called. It stores single condition UnivariateConditioned instances,
    for all conditions provided in cset argument, in addition to 'master' (no
    condition).

    Parameters
    ----------
    exp : :class:`Experiment` instance
        Defines which experiment, and how to parse it (with its filter set
        fset attribute)
    obs : :class:`Observable` instance
    eval_times : 1d ndarray
        array of times where statistics are evaluated
    cset : list of :class:`FilterSet` instances (default [])
        default: only 'master' is computed
    region : :class:`Region` instance or str (default 'ALL')
        with tmin, tmax, name attributes; first two are used to set boundaries
        for absolute time data (accept only acquisitions within those bounds),
        'name' is used for storage

    Attributes
    ----------
    _items : dictionary
        keys are condition labels, values are
        class:`UnivariateConditioned` instances.

    Raises
    ------
    UnivariateInitError
        when a parameter is missing (prints out which parameter)
    """

    def __init__(self, exp, obs, eval_times=None, region='ALL', cset=[]):
        self.obs = obs
        self.exp = exp
        self.cset = cset
        if eval_times is None:
            self._read_eval_times()
            if self.eval_times is None:
                raise UnivariateInitError('missing eval_times')
        else:
            self.eval_times = eval_times
        self.cset = cset
        if isinstance(region, Region):
            self.region = region
        elif isinstance(region, str):
            regs = Regions(exp)
            self.region = regs.get(region)
        else:
            regs = Regions(exp)
            self.region = regs.get('ALL')
        # create as many nodes as there are conditions in cset
        self._items = {}
        self._condition_labels = []
        # Instantiate the master UnivariateConditioned (no condition)
        master = UnivariateConditioned(self, applied_filter=None)
        # store it as 'master'
        self._items['master'] = master
        self._condition_labels.append('master')
        # create for each provided condition in cset
        self._add_conditions(cset)
        return

    def _add_conditions(self, cset=[]):
        """Create UnivariateConditioned instances associated to items in cset

        Parameters
        ----------
        cset : list of :class:`FilterSet` instances
        """
        for cdt in cset:
            cdt_repr = '{}'.format(repr(cdt))
            self._condition_labels.append(cdt_repr)
            self._items[cdt_repr] = UnivariateConditioned(self, cdt)
        return

    def _read_eval_times(self):
        """Read master item to get time array

        This one can be applied after reading data so as to update eval_times
        """
        master = self.master
        logger.debug('Reading evaluation times from master')
        try:
            master.read_text()
            self.eval_times = master.time
            logger.debug('Evaluation times import successful')
        except (analysis.MissingFileError, analysis.MissingFolderError):
            msg = 'Nothing to read there. Think of computing instead'
            print(msg)
            self.eval_times = None
            logger.debug('Evaluation times import failure')
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

    def _get_obs_path(self, user_root=None, write=False):
        """Get observable path"""
        obs = self.obs
        exp = self.exp
        fset = self.exp.fset
        analysis_path = analysis.get_analysis_path(exp, user_abspath=user_root,
                                               write=write)
        res = analysis.get_filter_path(analysis_path, fset, write=write)
        index_filter, filter_path = res
        obs_path = analysis.get_observable_path(filter_path, obs, write=write)
        return obs_path

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

    def export_text(self, analysis_folder=None):
        """Export results to text files.

        Parameters
        ----------
        analysis_folder : str (default None)
            Path to the analysis folder; default is 'analysis' subfolder in
            experiment folder
        """
        # write each condition
        for key, val in self._items.items():
            val.write_text(analysis_folder)
        return

    def import_from_text(self, analysis_folder=None):
        """Set instance from text files.

        This reader needs that Observable, Experiment, Conditions are defined
        and must match text folders.

        Raises
        ------
        UnivariateIOError
            when any of the folder is not found, or indexify does not match

        """
        # read each condition
        logger.debug('Trying to import from text files...')
        try:
            for key, val in self._items.items():
                logger.debug('Reading {} ..'.format(key))
                val.read_text(analysis_folder)
                logger.debug('..successful')
        except (analysis.MissingFileError, analysis.MissingFolderError) as missing:
            logger.debug('..failure: file is missing')
            raise UnivariateIOError(missing)
        except analysis.MismatchFileError as err:
            logger.debug('..failure: filename mismatch at level {}'.format(err.level))
            if err.level == 'observable':
                raise ObservableNameError('Obs name taken by a different observable')
        # update eval_times
        self._read_eval_times()
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
        Column names are ['time_interval', 'counts', 'auto_correlation']

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
        # [('time_interval', 'f8'),('counts', 'i4'), ('autocorrelation', 'f8')]
        return

    @property
    def time(self):
        if self.array is not None:
            return self.array['time_interval']
        else:
            return None

    @property
    def counts(self):
        if self.array is not None:
            return self.array['counts']
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
        cdt_path = cdt_univ._get_path(write=False)
        # otherwise a text.MissingFolder will be raised here
        if self.array is None:
            print('Nothing to write')
            return
        ffmt = '%.8e'  # floating point numbers
        ifmt = '%d'  # integers
        item_path = os.path.join(cdt_path, self.basename + '.tsv')
        names = self.array.dtype.names
        header = '\t'.join(names)
        fmt = [ifmt if 'counts' in n_ else ffmt for n_ in names]
        np.savetxt(item_path, self.array, fmt=fmt,
                   delimiter='\t', comments='', header=header)
        return

    def read_text(self, path='.'):
        """Initialize object by reading text output."""
        # get observable path that should exist already
        cdt_univ = self.statunivariate.univariate[self.condition]
        cdt_path = cdt_univ._get_path(write=False)
        item_path = os.path.join(cdt_path, self.basename + '.tsv')
        if not os.path.exists(item_path):
            raise analysis.MissingFileError(item_path)
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
        counts = self.array['counts']
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
    options : :class:`tunacell.stats.utils.CompuParams` instance
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
    options : :class:`tunacell.stats.utils.CompuParams` instance
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
        self.exp = univariate.exp
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

    def _get_obs_path(self, user_root=None, write=False):
        return self.univariate._get_obs_path(user_root=user_root, write=write)

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
        except analysis.MissingFolderError:
            self.single.export_text(analysis_folder)
            for key, val in self._items.items():
                val.write_text(analysis_folder)
        # export dataframe as csv file
        if self.dataframe is not None:
            exp = self.univariate.exp
            fset = self.univariate.exp.fset
            analysis_path = analysis.get_analysis_path(exp, user_abspath=analysis_folder,
                                                   write=True)
            res = analysis.get_filter_path(analysis_path, fset, write=True)
            index_filter, filter_path = res
            basename = 'data_{}_{}'.format(self.label, self.obs.name)
            text_file = os.path.join(filter_path, basename + '.csv')
            self.dataframe.to_csv(text_file, index=False)
        return

    def import_from_text(self, analysis_folder=None):
        try:
            for key, val in self._items.items():
                val.read_text(analysis_folder)
            exp = self.univariate.exp
            fset = self.univariate.exp.fset
            analysis_path = analysis.get_analysis_path(exp, user_abspath=analysis_folder,
                                                   write=False)
            res = analysis.get_filter_path(analysis_path, fset, write=False)
            index_filter, filter_path = res
            basename = 'data_{}_{}'.format(self.label, self.obs.name)
            text_file = os.path.join(filter_path, basename + '.csv')
            df = pd.read_csv(text_file, index_col=False)
            # convert column dtypes
            for col_name in df.columns:
                dtype = _dtype_converter(col_name)
                if dtype is not None:
                    df[col_name] = df[col_name].astype(dtype)
            self.dataframe = df
        except (analysis.MissingFileError, analysis.MissingFolderError):
            raise StationaryUnivariateIOError
        return
