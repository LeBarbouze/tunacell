#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Core computational module: cross-correlations and auto-correlations.
"""
from __future__ import print_function

import numpy as np
import pandas as pd
from scipy.linalg import toeplitz, triu
from scipy.interpolate import interp1d


class NoValidTimes(ValueError):
    pass


# %% Single observable computation of the statistics of dynamics
def set_dynamics(iter_timeseries, single, eval_times):
    """Central function that perform computations.

    It first defines dictionaries where rounded times are keys and values are
    iteratively updated. Then dictionaries are read to produce arrays that fill
    Univariate single instance.

    Parameters
    ----------
    iter_timeseries : iterator over TimeSeries instances
    single : initialized Univariate instance
    eval_times : 1d ndarray
        times at which statistics are computed

    Note
    -----
    This function computes:
    * times : times spanned by entire data (vector, size N determined by data)
      index i, time t_i :math:`\{ t_i\}_{i=1 \dots n}`
    * count_one : sample count of one point functions (vector, size N)
      :math:`\{c_i\}_{i=1 \dots n}=\text{number of samples at time } t_i `
    * sample_average -- sample averages at associated times (vector, size N)
      :math:`\{\langle x(t_i) \rangle\}_{i=1 \dots n}`
    * count_two : sample count of two points functions (matrix, size N*N)
      :math:`\{c^(2)_{ij}\}_{i,\ j=1 \dots n}` represents the number of
      samples connecting :math:`t_i` to :math:`t_j`. It is stored as a
      matrix of dimensions :math:`n \times n`.
    * autocorr : autocorrelation function (matrix, size N*N)

    .. math::

       a^{(2)}_{ij} = \langle \left( x(t_i) - \langle x(t_i) \rangle) \times
                      \langle \left( x(t_j) - \langle x(t_j) \rangle)

    These arrays are then stored in each Univariate._items entry as
    UnivariateConditioned instances, one for each condition, plus one
    for the unconditioned data ('master').
    """
    if single.region.name == 'ALL':
        tmin = None
        tmax = None
    else:
        tmin = single.region.tmin
        tmax = single.region.tmax
    # compute statistics and register in dictionaries
    # 'all' refers to unconditioned statistics
    records = {}
    # add counters for each condition
    for condition_lab in single._condition_labels:
        rec = {}
        rec['ones'] = np.zeros(len(eval_times))
        rec['count_ones'] = np.zeros(len(eval_times), dtype=int)
        rec['twos'] = np.zeros((len(eval_times), len(eval_times)))
        rec['count_twos'] = np.zeros((len(eval_times), len(eval_times)),
                                     dtype=int)
        records[condition_lab] = rec
    for ts in iter_timeseries:
        # loop over registered conditions in TimeSeries instance
        for condition_lab in ts.selections.keys():
            coords = ts.use_condition(condition_label=condition_lab,
                                      sharp_tleft=tmin, sharp_tright=tmax)
            if len(coords.clear_x) == 0:
                continue
            x = coords.clear_x
            y = coords.clear_y
            update(x, y, eval_times, records[condition_lab])

    # read individual counters and build results as 1d and 2d arrays
    for ic, condition_lab in enumerate(single._condition_labels):
        # average values
        one = records[condition_lab]['ones']
        count_one = records[condition_lab]['count_ones']

        mean = np.zeros_like(one)
        ok = np.where(np.logical_not(count_one == 0))
        mean[ok] = one[ok]/count_one[ok]
        mean[count_one == 0] = np.nan

        # matrix of covariances
        two = records[condition_lab]['twos']
        count_two = records[condition_lab]['count_twos']

        outer = np.zeros_like(two)
        ok = np.where(np.logical_not(count_two == 0))
        outer[ok] = two[ok]/count_two[ok]
        outer[count_two == 0] = np.nan
        # correct for mean
        mouter = np.outer(mean, mean)
        autocov = outer - mouter

        array = np.zeros(len(eval_times), dtype=[('time', 'f8'),
                                                 ('counts', 'u8'),
                                                 ('average', 'f8'),
                                                 ('std_dev', 'f8')])

        array['time'] = eval_times
        array['counts'] = count_one
        array['average'] = mean
        array['std_dev'] = np.sqrt(np.diag(autocov))

        # retrieve UnivariateConditioned instance
        conditioned_single = single._items[condition_lab]
        # associate data
        conditioned_single.bind(array, count_two, autocov)
    return


def update(time_array, value_array, eval_times, rec):
    """Update counters one and two.

    Parameters
    ----------
    time_array : 1d ndarray
        array of times
    value_array : 1d ndarray
        array of values, same size as time_array
    eval_times : 1d ndarray
        times at which functions are evaluated by interpolation method
    rec: dict
        dictionary that stores arrays to count and save one-, and two-point
        estimates, and that are updated with 't, val' timeseries sample.
    """
    # impossible to interpolate if less than 2 points
    if len(time_array) == 0:
        return
    # clean NaNs from values
    ok = np.where(np.logical_not(np.isnan(value_array)))
    t = time_array[ok]
    val = value_array[ok]
    if len(t) == 0:
        return
    if len(t) == 1:
        tt = t[0]
        ok = np.where(eval_times == tt)
        arr = np.zeros(len(eval_times))
        arr[:] = np.nan
        arr[ok] = val[0]
    else:
        f = interp1d(t, val, kind='linear', assume_sorted=True,
                     bounds_error=False)
        arr = f(eval_times)
    # check whether it's not all NaNs
    if np.all(np.isnan(arr)):
        return
    # find where it's not NaNs
    ok = np.where(np.logical_not(np.isnan(arr)))
    one = rec['ones']
    one[ok] = one[ok] + arr[ok]
    count_one = rec['count_ones']
    count_one[ok] = count_one[ok] + 1
    # same for the outer product
    outer = np.outer(arr, arr)
    ok = np.where(np.logical_not(np.isnan(outer)))
    two = rec['twos']
    two[ok] = two[ok] + outer[ok]
    count_two = rec['count_twos']
    count_two[ok] = count_two[ok] + 1
    return


# %% Computation of the stationary autocorrelation function
def set_stationary_autocorrelation(iter_timeseries, univariate, stationary,
                                   tmin=None, tmax=None, adjust_mean='global',
                                   disjoint=True):
    """Computes autocorrelation for stationary processes.

    Using univariate and parsing iter_timeseries, it computes autocorrelation
    function in stationary hypothesis. Result is stored in stationary,
    a StationaryUnivariate instance, which lists
    StationaryUnivariateConditioned instances, for each condition carried in
    stationary and parsed TimeSeries instances.

    Parameters
    ----------
    iter_timeseries : iterator over TimeSeries instances
    univariate : Univariate instance
        used to access average values
    stationary : StationaryUnivariate instance
        the one to be filled
    tmin : float
        lower bound for stationarity time range
    tmax : float
        upper bound for stationarity time range
    adjust_mean : str {'global', 'local'}
        how to substract average values: globally, or locally;
        use globally when local statistics are not sufficient,
        use locally is local statistics are sufficient.
    
    Raises
    ------
    NoValidTimes
        when time bounds make empty list of evaluation times
    """
    recs = {}  # one record per condition (including 'master")

    # we need to extract some information from univariate object (time, mean)
    local_means = {}
    if tmin is None:
        tmin = univariate.region.tmin
    if tmax is None:
        tmax = univariate.region.tmax

    time = univariate.master.time  # same for every condition
    # time-lapse timing: reduce dimension already (save computational time)
    if univariate.obs.timing == 't':
        window = np.logical_and(time >= tmin, time <= tmax)
        kwargs = {}
    # cell-cycle like : use tmin tmax to bound cell time values
    else:
        # eval_times are either generation index, either real times
        # but 'time' filtering will be done at cell-cycle using sharp bounds
        window = np.array(len(time) * [True, ])
        kwargs = {'sharp_tleft': tmin, 'sharp_tright': tmax}
    eval_times = time[window]
    if len(eval_times) == 0:
        raise NoValidTimes
    time_intervals = eval_times - eval_times[0]
    for condition_lab in stationary._condition_labels:

        # get UnivariateConditioned instance
        cdt_univariate = univariate[condition_lab]
        ms = cdt_univariate.average[window]
        co = cdt_univariate.count_one[window]  # needed to compute agg mean
        # construct dict of means indexed over indexified times
        local_mean = np.zeros_like(eval_times)
        if adjust_mean == 'local':
            local_mean = ms
        elif adjust_mean == 'global':
            agg_mean = np.nansum(co * ms)/np.nansum(co)
            local_mean = agg_mean * np.ones(len(eval_times))
        local_means[condition_lab] = local_mean
        # initialize rec
        recs[condition_lab] = {'counts': np.zeros(len(time_intervals), dtype=int),
                               'first': np.zeros(len(time_intervals)),
                               'second': np.zeros(len(time_intervals))}

    # store values
    df = None
    dfs = []
    # loop through timeseries
    for ts in iter_timeseries:
        df = ts.to_dataframe(**kwargs)
        if univariate.obs.timing == 't':
            df = df[np.logical_and(df.time >= tmin, df.time <= tmax)]
        dfs.append(df)
        for condition_lab in ts.selections.keys():
            coords = ts.use_condition(condition_label=condition_lab, **kwargs)
            if len(coords.clear_x) == 0:
                continue
            t = coords.clear_x
            v = coords.clear_y
            # get only times within valid window
            if univariate.obs.timing == 't':
                boo = np.logical_and(t >= tmin, t <= tmax)
            else:
                boo = np.array(len(t) * [True, ])
            rec = recs[condition_lab]  # this is where results are recorded
            local_mean = local_means[condition_lab]  # local means
            # update correlation
            update_stationary(t[boo], v[boo], eval_times, local_mean, rec, disjoint)

    df = pd.concat(dfs, ignore_index=True)
    stationary.dataframe = df

    # update each StationaryUnivariateConditioned instance
    for condition_lab in stationary._condition_labels:

        rec = recs[condition_lab]

        array = np.zeros(len(time_intervals),
                         dtype=[('time_interval', 'f8'),
                                ('counts', 'u8'),
                                ('auto_correlation', 'f8'),
                                ('std_dev', 'f8')])
        array['time_interval'] = time_intervals
        counts = rec['counts']
        first = rec['first']
        second = rec['second']
        ok = np.where(counts > 0)
        where_nan = np.where(counts == 0)
        array['counts'] = counts
        array['auto_correlation'][ok] = first[ok]/counts[ok]
        array['auto_correlation'][where_nan] = np.nan
        array['std_dev'][ok] = np.sqrt(second[ok]/counts[ok] - (first[ok]/counts[ok])**2)
        array['std_dev'][where_nan] = np.nan

        cdt_stationary = stationary[condition_lab]

        cdt_stationary.array = array
    return


def update_stationary(time_array, value_array, eval_times, local_mean, record,
                      disjoint=True):
    """Update counts and correlation value for stationary autocorrelation

    Parameters
    ----------
    t : 1d ndarray
        time values of timeseries
    val : 1d ndarray
        values of timeseries, corresponding to times
    eval_times : 1d ndarray
        times at which functions are evaluated using interpolation method
        (times outside the scope of 't' get np.nan values)
    local_mean : 1d ndarray
        array of average values, computed over the whole time values
        'eval_times'
    record : dict
        dictionary that saves arrays of counts, first, and second moments,
        saved as 1d ndarray
    disjoint : bool {True, False}
        whether to take disjoint time segments to evaluate statistics. When it
        is set to True, disjoint segments provide independent samples (under
        the Markovian assumption for the 't, val' process)
    """
    if len(time_array) == 0:
        return
    ok = np.where(np.logical_not(np.isnan(value_array)))
    t = time_array[ok]
    val = value_array[ok]
    if len(t) == 0:
        return
    if len(t) == 1:
        tt = t[0]
        ok = np.where(eval_times == tt)
        arr = np.zeros(len(eval_times))
        arr[:] = np.nan
        arr[ok] = val[0] - local_mean[ok]
    else:
        f = interp1d(t, val, kind='linear', assume_sorted=True,
                     bounds_error=False)
        arr = f(eval_times) - local_mean
    # check that it's not all NaNs:
    if np.all(np.isnan(arr)):
        return
    # get first index for which non-nan value
    for offset in range(len(arr)):
        if not np.isnan(arr[offset]):
            break
    outer = np.outer(arr, arr)
    for index in range(len(arr)-1):
        d = np.diagonal(outer, offset=index)
        if disjoint:
            step = index + 1
        else:
            step = 1
        sl = slice(offset, len(d), step)
        ok = np.logical_not(np.isnan(d[sl]))
        record['counts'][index] += len(d[sl][ok])
        record['first'][index] += np.nansum(d[sl])
        record['second'][index] += np.nansum(d[sl]*d[sl])
    return


# %% CROSS-CORRELATIONS

# UPDATING CROSS-CORRELATION COMPUTATION
def set_crosscorrelation(iter_timeseries, row_univ, col_univ, two):
    """Central function that computes cross-correlation matrix.

    It first defines dictionaries where rounded times are keys and values are
    iteratively updated. Then dictionaries are read to produce arrays that fill
    Univariate dynamics instance.

    Parameters
    ----------
    iter_timeseries : iterator over couple of TimeSeries instances
        see utils.iter_timeseries_2
    row_univ : :class:`Univariate` instance
        univariate instance corresponding to the first observable :math:`x`
    col_univ : :class:`Univariate` instance
        univariate instance corresponding to the second observable :math:`y`
    two : TwoObservable instance

    Note
    -----
    This function computes the cross-correlation matrix of observables (x, y):

    .. math::

       c^{(2)}_{ij} = \langle \left( x(t_i) - langle x(t_i) \rangle) \times
                      \langle \left( y(t_j) - \langle y(t_j) \rangle)

    These arrays are then stored in each TwoObservable._items entry as
    TwoObservableConditioned instances, one for each condition, plus one
    for the unconditioned data ('master').

    """
    # set condition list that match between both single instances
    cdt_labs = []
    for cdt_lab in row_univ._condition_labels:
        if cdt_lab in col_univ._condition_labels:
            cdt_labs.append(cdt_lab)

    # compute statistics and register in dictionaries
    # 'all' refers to unconditioned statistics
    records = {}
    means = {}

    row_eval_times = row_univ.eval_times
    col_eval_times = col_univ.eval_times
    # add counters for each common condition and define means
    for condition_lab in cdt_labs:
        # couple of dict indexified time : average value
        rec = {}
        means[condition_lab] = {'row': row_univ[condition_lab].average,
                                'col': col_univ[condition_lab].average}
        rec['counts'] = np.zeros((len(row_eval_times), len(col_eval_times)),
                                 dtype=int)
        rec['cross'] = np.zeros((len(row_eval_times), len(col_eval_times)))
        rec['square'] = np.zeros((len(row_eval_times), len(col_eval_times)))
        records[condition_lab] = rec

    for row_ts, col_ts in iter_timeseries:
        # loop over registered conditions in TimeSeries instance
        for condition_lab in cdt_labs:
            row_coords = row_ts.use_condition(condition_label=condition_lab,
                                              sharp_tleft=row_univ.region.tmin,
                                              sharp_tright=row_univ.region.tmax)
            col_coords = col_ts.use_condition(condition_label=condition_lab,
                                              sharp_tleft=col_univ.region.tmin,
                                              sharp_tright=col_univ.region.tmax)
            rec = records[condition_lab]
            row_mean = means[condition_lab]['row']
            col_mean = means[condition_lab]['col']
            update_2(row_coords, col_coords, row_eval_times, col_eval_times,
                     row_mean, col_mean, rec)

    # read individual counters and build results as 2d arrays
    for ic, condition_lab in enumerate(cdt_labs):
        # cross-correlation analysis
        rec = records[condition_lab]
        counts = rec['counts']
        cross = rec['cross']
        square = rec['square']
        corr = np.zeros(np.shape(cross))
        std = np.zeros(np.shape(cross))
        corr[counts == 0] = np.nan
        std[counts == 0] = np.nan
        corr[counts != 0] = cross[counts != 0]/counts[counts != 0]
        std[counts != 0] = np.sqrt(square[counts != 0]/counts[counts != 0] -
                                   (corr[counts != 0])*(corr[counts != 0]))
        condition_two = two._items[condition_lab]
        condition_two.counts = counts
        condition_two.cross = corr
        condition_two.std_dev = std
    return


def update_2(row_coords, col_coords, row_eval_times, col_eval_times,
             row_mean, col_mean, record):
    """Update counters one and two.

    Parameters
    ----------
    row_timeseries : Coordinates instance
        corresponding to first ('row') observable
    col_timeseries : Coordinates instance
        corresponding to second ('column') observable
    row_eval_times : 1d ndarray
        times at which row_timeseries is evaluated
    col_eval_times : 1d ndarray
        times at which col_timeseries is evaluated
    row_mean : 1d ndarray
        array of average values for row obs (same length as row_eval_times)
    col_mean : 1d ndarray
        array of average values for col obs (same length as col_eval_times)
    record : dict
        keys are: 'counts', 'cross', 'square'
        values are 2d ndarrays shape (len(row_eval_times), len(col_eval_times))

    Returns
    -------
    Update record entries
    """
    if len(row_coords.clear_x) == 0 or len(col_coords.clear_x) == 0:
        return
    # length 1 : take only the value if in eval_times
    if len(row_coords.clear_y) == 1:
        t, val = row_coords.clear_x[0], row_coords.clear_y[0]
        ok = np.where(row_eval_times == t)
        row_arr = np.zeros(len(row_eval_times))
        row_arr[:] = np.nan  # all NaNs but one
        row_arr[ok] = val - row_mean[ok]
    else:
        row_f = interp1d(row_coords.clear_x, row_coords.clear_y, kind='linear',
                         assume_sorted=True, bounds_error=False)
        row_arr = row_f(row_eval_times) - row_mean
    # if all NaNs, nothing to do
    if np.all(np.isnan(row_arr)):
        return
    if len(col_coords.clear_x) == 1:
        t, val = col_coords.clear_x[0], col_coords.clear_y[0]
        ok = np.where(col_eval_times == t)
        col_arr = np.zeros(len(col_eval_times))
        col_arr[:] = np.nan  # all NaNs but one
        col_arr[ok] = val - col_mean[ok]
    else:
        col_f = interp1d(col_coords.clear_x, col_coords.clear_y, kind='linear',
                         assume_sorted=True, bounds_error=False)
        col_arr = col_f(col_eval_times) - col_mean
    # if all NaNs, nothing to do
    if np.all(np.isnan(col_arr)):
        return
    out = np.outer(row_arr, col_arr)
    ok = np.where(np.logical_not(np.isnan(out)))
    record['counts'][ok] += 1
    record['cross'][ok] += out[ok]
    record['square'][ok] += out[ok] * out[ok]
    return


def set_stationary_crosscorrelation(iter_timeseries,
                                    row_univariate, col_univariate, stationary,
                                    tmin=None, tmax=None,
                                    adjust_mean='global',
                                    disjoint=True):
    """Edit stationary cross-correlation values from univariates and timeseries

    Parameters
    ----------
    iter_timeseries : iterator over couples of TimeSeries instances
        for each couple, first element corresponds to 'row', second to 'col'
    row_univariate : :class:`Univariate` instance
        corresponding to row observable, first items in iter_timeseries
    col_univariate : :class:`Univariate` instance
        corresponding to column observable, second items in iter_timeseries
    tmin : float (default None)
        lower bound (inclusive) for accepting time values
    tmax : float (default None)
        upper bound (inclusive) for accepting time values
    adjust_mean : str {'global', 'local'}
        option for computing mean value, substracted to all values
    disjoint : bool {True, False}
        whether to consider only disjoint time segments for sampling
    """
    # set condition list that match between both single instances
    cdt_labs = []
    for cdt_lab in row_univariate._condition_labels:
        if cdt_lab in col_univariate._condition_labels:
            cdt_labs.append(cdt_lab)

    if tmin is None:
        tmin = -np.infty
    if tmax is None:
        tmax = np.infty

    booarr = row_univariate.eval_times == col_univariate.eval_times
    if not np.all(booarr):
        raise ValueError('Evaluation times do not match!')
    time = row_univariate.eval_times  # same for every condition
    # time-lapse timing: reduce dimension already (save computational time)
    if row_univariate.obs.timing == 't':
        window = np.logical_and(time >= tmin, time <= tmax)
        kwargs = {}
    # cell-cycle like : use tmin tmax to bound cell time values
    else:
        # eval_times are either generation index, either real times
        # but 'time' filtering will be done at cell-cycle using sharp bounds
        window = np.array(len(time) * [True, ])
        kwargs = {'sharp_tleft': tmin, 'sharp_tright': tmax}
    eval_times = time[window]

    bwd = eval_times - eval_times[-1]
    fwd = eval_times - eval_times[0]
    time_intervals = np.concatenate([bwd[:-1], fwd])
    # print(eval_times)
    # print(time_intervals)

    means = {}
    recs = {}

    for cdt_lab in cdt_labs:
        means[cdt_lab] = {}
        for index, univariate in [('row', row_univariate),
                                  ('col', col_univariate)]:

            # get UnivariateConditioned instance
            cdt_univariate = univariate[cdt_lab]
            ms = cdt_univariate.average[window]
            co = cdt_univariate.count_one[window]  # needed to compute agg mean
            # construct dict of means indexed over indexified times
            local_mean = np.zeros_like(eval_times)
            if adjust_mean == 'local':
                local_mean = ms[window]
            elif adjust_mean == 'global':
                agg_mean = np.nansum(co * ms)/np.nansum(co)
                local_mean = agg_mean * np.ones(len(eval_times))
            means[cdt_lab][index] = local_mean
            # print(local_mean)

        # initialize rec
        recs[cdt_lab] = {'counts': np.zeros(len(time_intervals), dtype=int),
                         'first': np.zeros(len(time_intervals)),
                         'second': np.zeros(len(time_intervals))}

    # store values
    df = None
    dfs = []

    # loop through timeseries
    for row_ts, col_ts in iter_timeseries:
        # convert to dataframe first timeseries
        row_df = row_ts.to_dataframe(**kwargs)
        col_df = col_ts.to_dataframe(**kwargs)
        df = pd.merge(row_df, col_df, how='outer')
        # clean time values out of bounds
        if row_univariate.obs.timing == 't':
            df = df[np.logical_and(df.time >= tmin, df.time < tmax)]
        dfs.append(df)

        for condition_lab in cdt_labs:
            row_coords = row_ts.use_condition(condition_label=condition_lab,
                                              sharp_tleft=tmin,
                                              sharp_tright=tmax)
            col_coords = col_ts.use_condition(condition_label=condition_lab,
                                              sharp_tleft=tmin,
                                              sharp_tright=tmax)
            rec = recs[condition_lab]  # this is where results are recorded
            row_mean = means[condition_lab]['row']
            col_mean = means[condition_lab]['col']
            # update correlation
            update_stationary_cross(row_coords, col_coords, eval_times,
                                    row_mean, col_mean, rec,
                                    disjoint=disjoint)
    df = pd.concat(dfs, ignore_index=True)
    stationary.dataframe = df

    # update each StationaryUnivariateConditioned instance
    for condition_lab in cdt_labs:

        rec = recs[condition_lab]

        array = np.zeros(len(time_intervals), dtype=[('time_interval', 'f8'),
                                                     ('counts', 'u8'),
                                                     ('cross_correlation', 'f8'),
                                                     ('std_dev', 'f8')])
        array['time_interval'] = time_intervals
        counts = rec['counts']
        first = rec['first']
        second = rec['second']
        ok = np.where(counts > 0)
        where_nan = np.where(counts == 0)
        array['counts'] = counts
        array['cross_correlation'][ok] = first[ok]/counts[ok]
        array['cross_correlation'][where_nan] = np.nan
        array['std_dev'][ok] = np.sqrt(second[ok]/counts[ok] - (first[ok]/counts[ok])**2)
        array['std_dev'][where_nan] = np.nan

        cdt_stationary = stationary[condition_lab]

        cdt_stationary.array = array
    return


def update_stationary_cross(row_coords, col_coords, eval_times,
                            row_mean, col_mean, record, disjoint=True):
    """Update counts and correlation value for stationary cross-correlation

    Parameters
    ----------

    """
    if len(row_coords.clear_x) == 0 or len(col_coords.clear_x) == 0:
        return
    # length 1 : take only the value if in eval_times
    if len(row_coords.clear_y) == 1:
        t, val = row_coords.clear_x[0], row_coords.clear_y[0]
        ok = np.where(eval_times == t)
        row_arr = np.zeros(len(eval_times))
        row_arr[:] = np.nan  # all NaNs but one
        row_arr[ok] = val - row_mean[ok]
    else:
        row_f = interp1d(row_coords.clear_x, row_coords.clear_y, kind='linear',
                         assume_sorted=True, bounds_error=False)
        row_arr = row_f(eval_times) - row_mean
    # if all NaNs, nothing to do
    if np.all(np.isnan(row_arr)):
        return
    if len(col_coords.clear_x) == 1:
        t, val = col_coords.clear_x[0], col_coords.clear_y[0]
        ok = np.where(eval_times == t)
        col_arr = np.zeros(len(eval_times))
        col_arr[:] = np.nan  # all NaNs but one
        col_arr[ok] = val - col_mean[ok]
    else:
        col_f = interp1d(col_coords.clear_x, col_coords.clear_y, kind='linear',
                         assume_sorted=True, bounds_error=False)
        col_arr = col_f(eval_times) - col_mean
    # if all NaNs, nothing to do
    if np.all(np.isnan(col_arr)):
        return
    # get row first index for which non-nan value
    for row_offset in range(len(row_arr)):
        if not np.isnan(row_arr[row_offset]):
            break
    # get column first index for which non-nan value
    for col_offset in range(len(col_arr)):
        if not np.isnan(col_arr[col_offset]):
            break
    # find non-nan square matrix
    offset = np.amax([row_offset, col_offset])
    outer = np.outer(row_arr, col_arr)
    # eval_times : 0, 1, ..., n
    # time_intervals : -n, -n+1, ..., -1, 0, 1, ..., n
    n = len(eval_times) - 1
    # record for each time interval
    for index in np.arange(-n, n + 1):
        d = np.diagonal(outer, offset=index)
        if disjoint:
            step = np.abs(index) + 1
        else:
            step = 1
        sl = slice(offset, len(d), step)
        ok = np.logical_not(np.isnan(d[sl]))
        record['counts'][n + index] += len(d[sl][ok])
        record['first'][n + index] += np.nansum(d[sl])
        record['second'][n + index] += np.nansum(d[sl]*d[sl])
    return


def get_stat_from_dynamics(singlecdt, tmin=None, tmax=None):
    """Computes stationary autocorrelation vector from autocorr matrix.

    This function uses the autocorrelation matrix, that stores two-point
    autocorrelation functions between the various acquisition time-points,
    to compute the autocorrelation function in case of stationary hypothesis.
    Computation is fast but result is mostly unreliable (depends very much
    on the accuracy of the dynamical autocorrelation estimates, which is
    usually quite low). Use set_stationary_autocorrelation instead.

    Parameters
    ----------
    singlecdt : UnivariateConditioned instance
    tmin : float (default None)
    tmax : float (default None)

    Returns
    -------
    dts, cts, res
    dts : array of floats
        time intervals
    cts : array of ints
        sample counts
    res : array of floats
        autocorrelation values

    Note
    ----
    The estimate of the autocorrelation function using this procedure gives
    very poor accuracy estimates, and should be used only for quick inspection
    when a Univariate has been created and computed.
    For a better autocorrelation function estimate, it is necessary to parse
    samples another time, using only the sample average estimated in
    Univariate conditioned instances.
    """
    times = singlecdt.time
    autocorr = singlecdt.autocorr
    cts = singlecdt.count_two
    # Resize matrices depending on time limits
    indexlow, indexup = 0, None
    if tmin is not None:
        while indexlow < len(times) and times[indexlow] < tmin:
            indexlow += 1
    if tmax is not None:
        indexup = indexlow
        while indexup < len(times) and times[indexup] < tmax:
            indexup += 1
    sl = slice(indexlow, indexup)
    times = times[sl]
    autocorr = autocorr[sl, sl]
    cts = cts[sl, sl]
    # how many time-points
    nframes = len(times)
    all_counts = np.zeros(nframes, dtype=np.int)
    res = np.zeros(nframes, dtype=np.float)
    dts = np.zeros(nframes, dtype=np.float)
    col = np.zeros(nframes, dtype=np.int16)
    col[-1] = 1  # initialisation
    for k in range(nframes):
        col[k] = 1
        col[k-1] -= 1
        forward = triu(toeplitz(col))
        all_counts[k] = np.sum(forward * cts)
        res[k] = np.sum(forward * cts * autocorr)/all_counts[k]
        dts[k] = times[k] - times[0]
    return dts, all_counts, res
