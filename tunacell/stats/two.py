#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module defines objects storing results from cross-correlation analysis.
"""
from __future__ import print_function

import numpy as np
import os
import pandas as pd

from tunacell.io import analysis
from tunacell.stats.utils import _dtype_converter


class BivariateConditioned(object):
    """Stores dynamics conditioned statistics for a couple of observables.

    Parameters
    ----------
    obss : couple of :class:`Observable`  instances
    times : couple of 1d arrays
        each item refers to the array of times at which single obs conditioned
        statistics have been evaluated. The row and column indices of created
        matrices refers to the indices of first and second item respectively.
    applied_filter : :class:`FilterSet` instance
    """
    def __init__(self, bivariate, applied_filter=None):
        self.bivariate = bivariate
        self.applied_filter = applied_filter
        if applied_filter is not None:
            self.condition = repr(applied_filter)
        else:
            self.condition = 'master'
        # alias
        self.times = [uni.eval_times for uni in self.bivariate.univariates]
        self.counts = None  # 2-d array, sample counts at time t_i, t_j
        self.cross = None  # 2-d array, covariance matrix at (t_i, t_j)
        self.std_dev = None  # 2-d array, standard dev of covariances estimates
        return

    def as_dataframe(self):
        # unroll
        row_times = self.times[0]
        col_times = self.times[1]
        unroll_row_times = np.concatenate([np.array(len(col_times) * [rt, ]) for rt in row_times])
        unroll_col_times = np.concatenate([col_times for rt in row_times])
        counts = self.counts.flatten()
        cov = self.cross.flatten()
        std = self.std_dev.flatten()
        names = ['time-row', 'time-col', 'counts', 'covariance', 'std-cov']
        data = {'time-row': unroll_row_times,
                'time-col': unroll_col_times,
                'counts': counts,
                'covariance': cov,
                'std-cov' : std
                }
        return pd.DataFrame(data, columns=names)

    def _get_path(self, user_root=None, write=False):
        """Get condition path"""
        obs_path = self.bivariate._get_obs_path(user_root=user_root, write=write)
        res = analysis.get_condition_path(obs_path, self.applied_filter, write=write)
        index_condition, condition_path = res
        return condition_path

    def write_text(self, path=None):
        # 2 columns for times (already stored elsewhere, but just in case)
        cdt_path = self._get_path(user_root=path, write=True)
        item_path = os.path.join(cdt_path, 'times.tsv')
        with open(item_path, 'w') as f:
            f.write('row:')
            for t in self.times[0]:
                f.write('\t{:.2f}'.format(t))
            f.write('\n')
            f.write('column:')
            for t in self.times[1]:
                f.write('\t{:.2f}'.format(t))
            f.write('\n')
        # matrix for counts
        item_path = os.path.join(cdt_path, 'count_cross.tsv')
        np.savetxt(item_path, self.counts, fmt='%d', delimiter='\t')
        # matrix for cross-correlations
        item_path = os.path.join(cdt_path, 'cross.tsv')
        np.savetxt(item_path, self.cross, fmt='%.8e', delimiter='\t')
        item_path = os.path.join(cdt_path, 'std_dev.tsv')
        np.savetxt(item_path, self.std_dev, fmt='%.8e', delimiter='\t')
        return

    def read_text(self, path=None):
        cdt_path = self._get_path(user_root=path, write=False)
        item_path = os.path.join(cdt_path, 'times.tsv')
        if not os.path.exists(item_path):
            raise analysis.MissingFileError(item_path)
        times = []
        with open(item_path, 'r') as f:
            for line in f.readlines():
                times.append(list(map(float, line.rstrip().split('\t')[1:])))
        self.times = times
        # matrix for counts
        item_path = os.path.join(cdt_path, 'count_cross.tsv')
        if not os.path.exists(item_path):
            raise analysis.MissingFileError(item_path)
        arr = np.genfromtxt(item_path, dtype='i8', delimiter='\t')
        self.counts = arr
        # matrix for cross-correlations
        item_path = os.path.join(cdt_path, 'cross.tsv')
        if not os.path.exists(item_path):
            raise analysis.MissingFileError(item_path)
        arr = np.genfromtxt(item_path, dtype='f8', delimiter='\t')
        self.cross = arr
        # matrix for standard deviation of covariance estimates
        item_path = os.path.join(cdt_path, 'std_dev.tsv')
        if not os.path.exists(item_path):
            raise analysis.MissingFileError(item_path)
        arr = np.genfromtxt(item_path, dtype='f8', delimiter='\t')
        self.std_dev = arr
        pass

    def compute_stationary(self, indexify, tmin=None, tmax=None):
        """Computes stationary cross-correlation between tmin and tmax.

        DEPRECATED
        """
        # Narrow down valid times?
        time1 = self.single[str(self.obs[0])].time
        time2 = self.single[str(self.obs[1])].time

        rec = {}
        for index, t1 in enumerate(time1):
            if tmin is not None and t1 < tmin:
                continue
            if tmax is not None and t1 >= tmax:
                continue

            for jindex, t2 in enumerate(time2):
                if tmin is not None and t2 < tmin:
                    continue
                if tmax is not None and t2 >= tmax:
                    continue

                dt = t2 - t1
                sdt = indexify.indexify(dt)
                if sdt not in rec.keys():
                    rec[sdt] = [0, 0.]
                # Hum, I may have made a mistake, comment it:
                dcount = self.count_cross[index, jindex]
                rec[sdt][0] += dcount
                rec[sdt][1] += dcount * self.cross[index, jindex]

        dt_array = np.array(sorted(map(indexify.desindexify, rec.keys())))
        count_array = np.zeros(len(dt_array), dtype='u4')
        corr_array = np.zeros(len(dt_array), dtype='f8')

        for index, dt in enumerate(dt_array):
            sdt = indexify.indexify(dt)
            count_array[index] = rec[sdt][0]
            corr_array[index] = rec[sdt][1]/rec[sdt][0]

        array = np.zeros(len(dt_array), dtype=[('time_interval', 'f8'),
                                               ('count', 'u4'),
                                               ('cross-correlation', 'f8')])
        array['time_interval'] = dt_array
        array['count'] = count_array
        array['cross-correlation'] = corr_array

        return array


class BivariateError(Exception):
    pass


class BivariateIOError(IOError):
    pass


class Bivariate(object):
    """Stores dynamics statistics for a couple of observables.

    Parameters
    ----------
    row_univariate : :class:`Univariate`
        corresponds to row in cross-correlation matrices
    col_univariate : :class:`Univariate`
        corresponds to colum in cross-correlation matrices

    """
    def __init__(self, row_univariate, col_univariate):
        # check whether exp instances match
        s1, s2 = row_univariate, col_univariate
        if s1.exp.abspath != s2.exp.abspath:
            raise BivariateError('Experiments do not match')
        if repr(s1.exp.fset) != repr(s2.exp.fset):
            raise BivariateError('Filter sets do not match')
        self.univariates = (row_univariate, col_univariate)
        self.exp = s1.exp
        # build common conditions
        cset = []
        for cdt in s1.cset:
            if repr(cdt) in map(repr, s2.cset):
                cset.append(cdt)
        self.cset = cset
        self._items = {}
        self._condition_labels = ['master', ]
        # alias
        Bic = BivariateConditioned
        self._items['master'] = Bic(self, applied_filter=None)
        for cdt in cset:
            lab = repr(cdt)
            self._condition_labels.append(lab)
            self._items[lab] = Bic(self, applied_filter=cdt)
        return

    def _get_obs_path(self, user_root=None, write=False):
        """Get observable path"""
        obss = [univ.obs for univ in self.univariates]
        exp = self.exp
        fset = self.exp.fset
        analysis_path = analysis.get_analysis_path(exp, user_abspath=user_root,
                                               write=write)
        res = analysis.get_filter_path(analysis_path, fset, write=write)
        index_filter, filter_path = res
        obs_path = analysis.get_biobservable_path(filter_path, obss, write=write)
        return obs_path

    def export_text(self, analysis_folder=None):
        # write each condition
        for key, val in self._items.items():
            val.write_text(analysis_folder)
        return

    def import_from_text(self, analysis_folder=None):
        # read each condition
        try:
            for key, val in self._items.items():
                val.read_text(analysis_folder)
        except (analysis.MissingFileError, analysis.MissingFolderError) as missing:
            raise BivariateIOError(missing)
        return

    def __getitem__(self, key):
        return self._items[key]

    @property
    def master(self):
        """There's always a master (no condition)"""
        return self['master']


class StationaryBivariateConditioned(object):
    """Cross-correlation as a function of time difference, for a univariate cdtion.

    Parameters
    ----------
    obss : couple of :class:`Observable`  instances
    times : couple of 1d arrays
        each item refers to the array of times at which univariate obs conditioned
        statistics have been evaluated. The row and column indices of created
        matrices refers to the indices of first and second item respectively.
    applied_filter : :class:`FilterSet` instance
    tmin : float (default None)
    tmax : float (default None)
    adjust_mean : str {'global', 'local'}
    """

    def __init__(self, statbivariate, applied_filter=None, array=None):
        self.statbivariate = statbivariate
        self.basename = 'stationary_cross'
        # add region label
        self.basename += '_' + self.statbivariate.region.name
        # add computation options
        self.basename += '_' + self.statbivariate.options.as_string_code()
        self.applied_filter = applied_filter
        self.basename = 'stationary_bivariate'
        if applied_filter is not None:
            self.condition = repr(applied_filter)
        else:
            self.condition = 'master'
        self.array = array  # should be a 3 columns array
        return

    def as_dataframe(self):
        return pd.DataFrame(self.array)

    def _get_path(self, user_root=None, write=False):
        """Get condition path"""
        obs_path = self.statbivariate._get_obs_path(user_root=user_root, write=write)
        res = analysis.get_condition_path(obs_path, self.applied_filter, write=write)
        index_condition, condition_path = res
        return condition_path

    def write_text(self, path='.'):
        """Write array to file."""
        # get condition p
        cdt_path = self._get_path(user_root=path, write=True)
        if self.array is None:
            print('Nothing to write')
            return
        ffmt = '%.8e'  # floating point numbers
        ifmt = '%d'  # integers
        item_path = os.path.join(cdt_path, self.basename + '.tsv')
        names = self.array.dtype.names
        header = '\t'.join(names)
        fmt = [ifmt if 'count' in n_ else ffmt for n_ in names]
        np.savetxt(item_path, self.array, fmt=fmt,
                   delimiter='\t', comments='', header=header)
        return

    def read_text(self, path='.'):
        """Initialize object by reading text output."""
        cdt_path = self._get_path(user_root=path, write=False)
        item_path = os.path.join(cdt_path, self.basename + '.tsv')
        if not os.path.exists(item_path):
            raise analysis.MissingFileError(item_path)
        arr = np.genfromtxt(item_path, delimiter='\t', names=True)
        self.array = arr
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
    def crosscorr(self):
        if self.array is not None:
            return self.array['cross_correlation']
        else:
            return None


class StationaryBivariateIOError(IOError):
    pass


class StationaryBivariate(object):
    """Cross-correlation analysis at stationarity

    Parameters
    ----------
    row_univariate : :class:`Univariate` instance
    col_univariate : :class:`Univariate` instance
    region : object with tmin, tmax, name attributes
        determine where the process is (hypothetically) stationary
    options : :class:`CompuParams` instance
    """

    def __init__(self, row_univariate, col_univariate,
                 region=None, options=None):
        self.region = region
        self.options = options
        s1, s2 = row_univariate, col_univariate
        # obss = [univariate.obs for univariate in univariates]
        if s1.exp.abspath != s2.exp.abspath:
            raise BivariateError('Experiments do not match')
        if repr(s1.exp.fset) != repr(s2.exp.fset):
            raise BivariateError('Filter sets do not match')

        self.univariates = (row_univariate, col_univariate)
        self.label = self.region.name
        self.tmin = self.region.tmin
        self.tmax = self.region.tmax
        self.adjust_mean = self.options.adjust_mean
        self.disjoint = self.options.disjoint
        self.exp = s1.exp
        # build common conditions
        cset = []
        for cdt in s1.cset:
            if repr(cdt) in map(repr, s2.cset):
                cset.append(cdt)
        self.cset = cset
        self.dataframe = None  # to be updated with pandas.DataFrame object
        self._condition_labels = ['master', ]
        self._items = {}
        # alias
        SBic = StationaryBivariateConditioned
        self._items['master'] = SBic(self, applied_filter=None, array=None)

        for cdt in cset:
            self._condition_labels.append(repr(cdt))
            self._items[repr(cdt)] = SBic(self, applied_filter=cdt, array=None)
        return

    def _get_obs_path(self, user_root=None, write=False):
        """Get observable path"""
        obss = [univ.obs for univ in self.univariates]
        exp = self.exp
        fset = self.exp.fset
        analysis_path = analysis.get_analysis_path(exp, user_abspath=user_root,
                                               write=write)
        res = analysis.get_filter_path(analysis_path, fset, write=write)
        index_filter, filter_path = res
        obs_path = analysis.get_biobservable_path(filter_path, obss, write=write)
        return obs_path

    def export_text(self, analysis_folder=None):
        # write each condition
        try:
            for key, val in self._items.items():
                val.write_text(analysis_folder)
        # when not possible it means single object has not been exported yet
        except analysis.MissingFolderError:
            for uni in self.univariates:
                uni.export_text(analysis_folder)
            for key, val in self._items.items():
                val.write_text(analysis_folder)
        # export dataframe as csv file
        if self.dataframe is not None:
            exp = self.exp
            fset = self.exp.fset
            analysis_path = analysis.get_analysis_path(exp,
                                                   user_abspath=analysis_folder,
                                                   write=True)
            res = analysis.get_filter_path(analysis_path, fset, write=True)
            index_filter, filter_path = res
            o1, o2 = [uni.obs for uni in self.univariates]
            basename = 'data_{}_{}---{}'.format(self.label, o1.name, o2.name)
            text_file = os.path.join(filter_path, basename + '.csv')
            self.dataframe.to_csv(text_file, index=False)
        return

    def import_from_text(self, analysis_folder=None):
        try:
            for key, val in self._items.items():
                val.read_text(analysis_folder)
            exp = self.exp
            fset = self.exp.fset
            analysis_path = analysis.get_analysis_path(exp, user_abspath=analysis_folder,
                                                   write=False)
            res = analysis.get_filter_path(analysis_path, fset, write=False)
            index_filter, filter_path = res
            o1, o2 = [uni.obs for uni in self.univariates]
            basename = 'data_{}_{}---{}'.format(self.label, o1.name, o2.name)
            text_file = os.path.join(filter_path, basename + '.csv')
            if not os.path.exists(text_file):
                raise analysis.MissingFileError
            df = pd.read_csv(text_file, index_col=False)
            # convert column dtypes
            for col_name in df.columns:
                dtype = _dtype_converter(col_name)
                if dtype is not None:
                    df[col_name] = df[col_name].astype(dtype)
            self.dataframe = df
        except (analysis.MissingFileError, analysis.MissingFolderError):
            raise StationaryBivariateIOError
        return

    def __getitem__(self, key):
        return self._items[key]

    @property
    def master(self):
        """There's always a master (no condition)"""
        return self['master']
