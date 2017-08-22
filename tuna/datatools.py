#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
This module implements general data operations.
"""
from __future__ import print_function  # start to adapt to Python 3

import numpy as np
from scipy.interpolate import interp1d
from numpy.lib.recfunctions import append_fields
import warnings


class MissingLabel(Exception):
    pass


class LocalFitError(Exception):
    pass


class Coordinates(object):
    """Clear timeseries from NaNs and other functionalities"""

    def __init__(self, x, y):
        if np.any(np.isnan(x)):
            msg = 'NaN value(s) detected in x-array: please remove beforehand'
            raise ValueError(msg)
        if len(x) != len(y):
            msg = 'Arrays of different lengths! Check x and y input'
            raise ValueError(msg)
        self._x = x
        self._y = y
        self.valid = np.logical_not(np.isnan(y))
        return

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def clear_x(self):
        if len(self._y) > 0:
            return self._x[self.valid]
        else:
            return np.array([], dtype=float)

    @property
    def clear_y(self):
        if len(self._y) > 0:
            return self._y[self.valid]
        else:
            return np.array([], dtype=float)


# %% NEW LOCAL FIT ESTIMATE USING ARRAYS

def compute_rates(x, y, x_break=None,
                  anterior_x=[], anterior_y=[],
                  scale='log',
                  time_window=15., dt=5.,
                  join_points=3,
                  testing=False):
    """Computes rates of array y against x by local fits over shifting window.

    Results are evaluated at coordinate array x by linear interpolation, so
    that fitted and rate arrays are of same size of x and y.
    When possible, it uses anterior values (in parent cell in tunacell context)
    to extend the result range. When a timeseries is expected to cover multiple
    cell cycles, with continuity hypothesis at the level of rates, it allows to
    use parent cell data to extend the range over which the dferivative can be
    evaluated.

    Parameters
    ----------
    x : 1d ndarray
        co-ordinate array (usually array of times for timeseries)
    y : 1d ndarray
        ordinate (array of values of same length as x array)
    x_break : float
        value of co-ordinate at which continuity joining is performed: in
        tunacell context, this is the time of birth for present cell, and it
        corresponds to division time of its parent cell
    anterior_x : 1d ndarray
        extra-bound anterior co-ordinate array: in tunacell context, this is
        the time array extracted from parent cell
    anterior_y : 1d ndarray
        extra-bound anterior ordinate array: in tunacell context, this is the
        array of values extracted from parent cell
    scale : str {'linear', 'log'}
        expected scale of y versus x. For exponential growth, use 'log' scale.
    time_window : float
        size of time window over which local fit is performed
    dt : float
        acquisition period of time array
    join_points : int (default 3)
        minimal number of points used when performing local fits to make
        continuity between anterior and present timeseries
    testing : bool {False, True}
        verbose output for testing

    Returns
    -------
    rates : 1d ndarray
        array of rates estimated at x array (uses interpolation)
    fits : 1d ndarray
        array of fitted values at x array
    used_x : 1d ndarray
        concatenation with derivative continuity hypothesis of anterior_x and x
    used_y : 1d ndarray
        concatenation with derivative continuity hypothesis of anterior_y and y
    """
    # check array lengths
    coords = Coordinates(x, y)
    anterior_coords = Coordinates(anterior_x, anterior_y)

    # define x_break by convention
    if x_break is None:
        x_break = coords.x[0] - dt/2.  # artificial birth time (CONVENTION)

    #  operators and their inverse
    if scale == 'log':
        y_operator = np.log
        y_inv_operator = np.exp
    elif scale == 'linear':
        def y_operator(vals):
            return vals

        def y_inv_operator(vals):
            return vals

    # remove NaN values
    clean_x = coords.clear_x
    clean_y = coords.clear_y
    clean_ax = anterior_coords.clear_x
    clean_ay = anterior_coords.clear_y

    # find period (can be different from dt when multiple acquisition periods)
    if len(clean_x) == 0:
        return [], [], [], []
    # too small x range to fit anything
    elif clean_x[-1] - clean_x[0] < time_window:
        only_nans = len(x) * [np.nan, ]
        return only_nans, only_nans, x, y

    # clean_x is necesarily of length >=2
    period = np.amin(clean_x[1:] - clean_x[:-1])

    # number of points per time window
    n_points = int(np.round(time_window/period, decimals=0))
    if testing:
        print('Local fits will be performed over {} points'.format(n_points))
    if n_points < 2:
        msg = ('Trying to perform linear fit over less than 2 points. '
               'Please use a larger time_window parameter.')
        raise LocalFitError(msg)
    elif n_points == 2:
        msg = ('Performing linear fit over 2 points: '
               'for rate computation experimental errors are not smoothed')
        warnings.warn(msg)

    # define the number of points for estimates at cell birth, parent division
    if n_points >= join_points:
        n_joints = n_points
    else:
        n_joints = join_points
    if testing:
        print('For values at division/birth:')
        print('  try to fit over {} points (local fits)'.format(n_points))
        print('  but allows to reduce up to {} points'.format(
                join_points))
        print()

    # auxiliary variables
    op_y_break = None  # estimate of y value at joining (value at birth)
    op_ay_break = None  # estimate of anterior y value at joining (at division)

    op_y = y_operator(clean_y)

    if len(clean_x) >= join_points:
        # fit to at least join_points, more if possible
        r, i = np.polyfit(clean_x[:n_joints], op_y[:n_joints], 1)
        op_y_break = i + r * x_break

    if testing:
        msg = ('Data to fit:\n'
               'x : {}'.format(clean_x) + '\n'
               'y : {}'.format(clean_y))
        print(msg)

    # try to use anterior values : compute break offset
    trans_op_ay = []  # translated, operated anterior values; default: empty
    if len(clean_ax) > 0:
        op_ay = y_operator(clean_ay)
        # 2 checks:
        #   1. there at enough points to get the final value estimate
        cdt1 = len(clean_ax) >= join_points
        #   2. initial value is determined
        cdt2 = op_y_break is not None
        if cdt1 and cdt2:
            r, i = np.polyfit(clean_ax[-n_joints:],
                              op_ay[-n_joints:], 1)
            op_ay_break = i + r * x_break
            if testing:
                msg = ('Extrapolated anterior value at break:\n '
                       'break time, value: '
                       '{}, {}'.format(x_break, y_inv_operator(op_ay_break)))
                print(msg)
                print()
            # adjust values by translating
            trans_op_ay = op_ay - op_ay_break + op_y_break
    if testing:
        msg = ('Translated, operated anterior values are:\n'
               '{}'.format(trans_op_ay))
        print(msg)

    # concatenate operated values
    all_x = np.concatenate([clean_ax, clean_x])
    all_op_y = np.concatenate([trans_op_ay, op_y])
    all_y = y_inv_operator(all_op_y)

    fit_x = np.zeros_like(all_x)
    fit_op_y = np.zeros_like(all_x)  # fitted values
    rate_op_y = np.zeros_like(all_x)   # rates of local fits

    # sliding window
    for index, t in enumerate(all_x):
        t_start = t - period/2.  # convention
        t_stop = t_start + time_window
        # time of evaluation
        time_eval = (t_start + t_stop)/2.
        fit_x[index] = time_eval
        # check that at least one time point is larger than break point
        if t_stop <= x_break:
            fit_op_y[index] = np.nan
            rate_op_y[index] = np.nan
            continue
        # reduce to points to fit
        lower = all_x > t_start
        upper = all_x <= t_stop
        boo = np.logical_and(lower, upper)
        if testing:
            print('+ Time window:', end=' ')
            print('{} < t < {} '.format(t_start, t_stop), end=' ')
            print('({} points)'.format(len(all_x[boo])))
        # check if n_points > 3
        if len(all_x[boo]) < n_points:
            if testing:
                print('Not enough points on this time window')
                print('({} instead of {})'.format(len(all_x[boo]), n_points))
                print('Going to next time window')
                print()
            # not enough point in this window: insert NaN
            rate_op_y[index] = np.nan
            fit_op_y[index] = np.nan
        else:
            if testing:
                print('local fit over {} points'.format(len(all_x[boo])))
            rate, intercept = np.polyfit(all_x[boo], all_op_y[boo], 1)
            fit_op_y[index] = rate * time_eval + intercept
            rate_op_y[index] = rate

    # new values: linear interpolation
    fit_op_coords = Coordinates(fit_x, fit_op_y)
    f = interp1d(fit_op_coords.clear_x, fit_op_coords.clear_y, kind='linear',
                 assume_sorted=True, bounds_error=False)
    out_y = np.array(len(coords.x) * [np.nan, ])  # initialize to NaNs
    # valid data points: interpolation at initial x coordinates
    out_y[coords.valid] = y_inv_operator(f(coords.clear_x))

    # rates : linear interpolation
    rate_coords = Coordinates(fit_x, rate_op_y)
    f = interp1d(rate_coords.clear_x, rate_coords.clear_y, kind='linear',
                 assume_sorted=True, bounds_error=False)
    out_rate = np.array(len(coords.x) * [np.nan, ])
    # valid data points: interpolation at initial x coordinates
    out_rate[coords.valid] = f(coords.clear_x)

    return out_rate, out_y, all_x, all_y


# %% build timeseries over cell instances

def local_rate(cell, yaxis='length', yscale='log',
               time_window=15., dt=5.,
               join_points=3,
               testing=False):
    """Performs derivative estimate of observable yaxis over variable 'time'.

    Shifting time window estimate of 'time' -> d(yaxis)/d('time') (axis)
    using cell and parent cell values when possible.

    Use this function when you expect that the derivative (of log) of yaxis is
    continuous at cell divisions.

    Parameters
    ----------
    cell : Cell instance
       (should have a .data attribute, being a Numpy structured array)
    yaxis : str (default 'length')
       member of cell.data.dtypes.names
       (measurements of this variable are expected to be noisy, hence the local
       fits)
    yscale : str (default 'log'')
       how data is treated to join parent end/cell start, and for rates
       (use 'log' for a quantity that you expect to be exponentially growing,
       otherwise 'linear' is a safe call)
    time_window : float
       time window over which data points are collected for fitting procedure
    dt : float
       expected 'time' interval between two consecutive data points in array
    join_points : int (default 3)
       minimal number of data points taken to estimate cell cycle boundary
       values (birth value of cell, division value of parent cell)
    testing : bool (default False)
       set to True to have a verbose description printed on standard output
       while function is called.

    Returns
    -------
    output, adjusted_yaxis_timeseries

    output : NumPy structured array
       columns: time, cellID, rate_<yaxis>, fit_<yaxis>
       rate_<yaxis> is the derivative of (optionally log of) yaxis, estimated
       by a the linear coefficient of the fit over time_window
       fit_<yaxis>: fit over the time window of <yaxis> values

    adjusted_yaxis : NumPy structured array
        columns: time, cellID, <yaxis> label
        corresponds to the concatenation of parent and daughter values, with
        parent values adjusted to ensure 'continuity' at division

    Notes
    -----
    Try to put time_window as a multiple of dt, and avoid values of the from
    (2q+1)*dt/2, since it may create undesired boundaries effect.
    """
    # derivative of values, or of the log(values)
    if yscale == 'log':
        y_operator = np.log
        y_inv_operator = np.exp
    elif yscale == 'linear':
        y_operator = lambda x: x
        y_inv_operator = lambda x: x

    output = []
    yaxis_timeseries = []

    if testing:
        print('Starting to compute local growth rates of ' + yaxis)
        print('for cell: {}'.format(cell))

    epsilon = dt/10.  # auxiliary delta time

    # number of points per time window
    n_points = int(np.round(time_window/dt, decimals=0))
    if testing:
        print('Local fits will be performed over {} points'.format(n_points))
        
    if n_points < 2:
        msg = ('Trying to perform linear fit over less than 2 points. '
               'Please use a larger time_window parameter.')
        raise LocalFitError(msg)
    elif n_points == 2:
        msg = ('Performing linear fit over 2 points: '
               'for rate computation experimental errors are not smoothed')
        warnings.warn(msg)
    value_0 = None  # estimate of cell log initial value
    T0 = None  # start of time for cell
    offset = None
    times = []
    values = []
    cids = []
    parent_times = []
    adjusted_parent_values = []
    pids = []
    id_type = 'u2'  # default

    # define the number of points for estimates at cell birth, parent division
    if n_points >= join_points:
        n_joints = n_points
    else:
        n_joints = join_points
    if testing:
        print('For values at division/birth:')
        print('  try to fit over {} points (local fits)'.format(n_points))
        print('  but allows to reduce up to {} points'.format(
                join_points))
        print()

    # find initial time and if possible initial value
    if cell.data is not None and len(cell.data) > 0:
        times = cell.data['time']
#        cids = cell.data['cellID']
#        id_type = cids.dtype
        values = y_operator(cell.data[yaxis])
        if cell.birth_time is not None:
            T0 = cell.birth_time
        else:
            T0 = times[0] - dt/2.  # artificial birth time
        if len(times) >= join_points:
            # fit to at least join_points, more if possible
            r, i = np.polyfit(times[:n_joints],
                              values[:n_joints], 1)
            value_0 = i + r * T0

        if testing:
            print('Okay we found some data for input cell:')
            print('times: {}'.format(times))
            print(yaxis + ': {}'.format(y_inv_operator(values)))
            if value_0 is not None:
                msg = ('We extrapolate this value for birth value'
                       'time, value: '
                       '{}, {}'.format(T0, y_inv_operator(value_0)))
                print(msg)
            print()

        # let's go: find parent values if they exist
        parent_value_1 = None  # estimate of parent end log value
        if cell.parent is not None and value_0 is not None:
            # use parent cell for first estimates
            if cell.parent.data is not None:
                PT1 = cell.parent.division_time  # parent time of division
                parent_times = cell.parent.data['time']
#                pids = cell.parent.data['cellID']
                parent_values = y_operator(cell.parent.data[yaxis])
                if testing:
                    print('Okay, we found some data for parent cell')
                    print('times: {}'.format(parent_times))
                    print(yaxis + ': {}'.format(y_inv_operator(parent_values)))
                # 2 checks:
                #   1. there at enough points to get the final value estimate
                cdt1 = len(parent_times) >= join_points
                #   2. check if parent division and cell division coincide
                cdt2 = np.abs(PT1 - T0) < epsilon
                #   3. intial value is determined
                cdt3 = value_0 is not None
                if cdt1 and cdt2 and cdt3:
                    r, i = np.polyfit(parent_times[-n_joints:],
                                      parent_values[-n_joints:], 1)
                    parent_value_1 = i + r * PT1
                    if testing:
                        msg = ('We extrapolate this value for division value'
                               'time, value: '
                               '{}, {}'.format(PT1,
                                               y_inv_operator(parent_value_1)))
                        print(msg)
                        print()
                    # adjust values by translating
                    offset = parent_value_1 - value_0
                    adjusted_parent_values = parent_values - offset
                    if testing:
                        msg = ('Values of parent cell are rescaled such that:'
                               'rescaled ' + yaxis + ': '
                               '{}'.format(adjusted_parent_values))
                        print(msg)
                # not enough point in parent cell to estimate parent_value_1
                # reset to empty lists
                else:
                    parent_times = []
#                    pids = []
                    adjusted_parent_values = []
                    # maybe not enough point in parent to find adjacent timepoint
#                    while np.amin(np.abs(parent_times - t_start)) > dt + epsilon:
#                        t_start += dt

        # concatenate arrays (may be done in structured arrays)
        ts = np.concatenate((parent_times, times))
        vals = np.concatenate((adjusted_parent_values, values))
#        ids = np.concatenate((pids, cids))
        if testing:
            print('size of ts {}, size of vals {}'.format(len(ts), len(vals)))

        # define adjusted timeseries as structured array
        label = 'adjusted_' + yaxis
        yaxis_timeseries = np.zeros(len(ts), dtype=[('time', 'f8'),
                                                    (label, 'f8'),
#                                                    ('cellID', id_type)
                                                    ]
                                    )
        yaxis_timeseries['time'] = ts
        yaxis_timeseries[label] = y_inv_operator(vals)
#        yaxis_timeseries['cellID'] = ids

#        yaxis_timeseries = zip(ts, y_inv_operator(vals))

        # start the earliest possible in parent cell such that
        # the time window intercepts at least one frame in cell
        t_start = T0 - time_window + dt
        t_stop = t_start + time_window
        # values within cell cycle

        if testing:
            print('Loop to compute local rates')

        # in principle, there are as many points as cell timepoints

        new_times = np.zeros_like(times)
#        new_ids = np.zeros(len(times), dtype=id_type)
        new_vals = np.zeros(len(new_times), dtype='f8')
        rates = np.zeros(len(new_times), dtype='f8')

        for index, t in enumerate(times):
            t_stop = t + dt/2.  # convention
            t_start = t_stop - time_window
            # time of evaluation
            time_eval = (t_start + t_stop)/2.
            new_times[index] = time_eval
            # reduce to points to fit
            lower = ts > t_start
            upper = ts <= t_stop
            boo = np.logical_and(lower, upper)
            if testing:
                print('+ Time window:', end=' ')
                print('{} < t < {} '.format(t_start, t_stop), end=' ')
                print('({} points)'.format(len(ts[boo])))
            # check if n_points > 3
            if len(ts[boo]) < n_points:
                if testing:
                    print('Not enough points on this time window')
                    print('({} instead of {})'.format(len(ts[boo]), n_points))
                    print('Going to next time window')
                    print()
                # not enough point in this window: insert NaN
                rates[index] = np.nan
                new_vals[index] = np.nan
            else:
                if testing:
                    print('local fit over {} points'.format(len(ts[boo])))
                rate, intercept = np.polyfit(ts[boo], vals[boo], 1)
                new_vals[index] = rate * time_eval + intercept
                rates[index] = rate
        # offset parent values
        boo = new_times < T0
        if offset is not None:
            new_vals[boo] = new_vals[boo] + offset
        # assign cellID column
#        if len(pids) > 0:
#            new_ids[boo] = pids[0]
#        new_ids[np.logical_not(boo)] = cids[0]

        # inverse operator for fitted values
        new_vals = y_inv_operator(new_vals)

        # build structured array for output
        output = np.zeros(len(new_times), dtype=[('time', 'f8'),
                                                 ('fit_' + yaxis, 'f8'),
                                                 ('rate_' + yaxis, 'f8'),
#                                                 ('cellID', id_type)
                                                 ])
        output['time'] = new_times
        output['fit_' + yaxis] = new_vals
        output['rate_' + yaxis] = rates
#        output['cellID'] = new_ids

        # clean NaNs
        boo = np.logical_not(np.isnan(rates))
        output = output[boo]

    return output, yaxis_timeseries


class ExtrapolationError(Exception):
    pass


# TODO : extend to rates generated observables using above function?
def extrapolate_endpoints(cell, timeseries=None, yaxis='length', scale='log',
                          end_point='birth', join_points=3):
    """Extrapolate values at cell birth and cell division.

    Parameters
    ----------
    cell : Cell instance
    timeseries : optional, sequence of (time, value) (default None)
        if provided, the extrapolation will use these values
        if not, the sequence of (time, value) is retrieved from cell.data
        using yaxis parameter
    yaxis : str
        raw observable to evaluate at end_point (e.g. 'length')
    scale : str, {'linear', 'log'}
        expected scale of observable
    end_point : str, {'birth', 'division'}
    join_points : int, number of points over which fit is done

    Returns
    -------
    value - value estimated at birth/division for <yaxis> observable

    Raises
    ------
    ExtrapolationError: when extrapolation fails due to too less points, or
                        when birth or division are not defined on Cell instance
    """
    npts = join_points
    if scale == 'log':
        y_operator = np.log
        y_inv_operator = np.exp
    elif scale == 'linear':
        y_operator = lambda x: x
        y_inv_operator = lambda x: x

    # arrays of time and observable
    if timeseries is not None:
        times, values = map(np.array, zip(*timeseries))
    else:
        times, values = map(np.array, zip(*cell.data[['time', yaxis]]))

    op_values = y_operator(values)
    if len(times) < join_points:
        msg = ('not enough points to fit end-point: '
               '{} instead of {}'.format(len(times), join_points))
        raise ExtrapolationError(msg)

    if end_point == 'birth' or end_point == 'b':
        if cell.birth_time is None:
            raise ExtrapolationError('cell with no birth time defined')
        rate, intercept = np.polyfit(times[:npts], op_values[:npts], 1)
        val = rate * cell.birth_time + intercept

    if end_point == 'division' or end_point == 'd':
        if cell.division_time is None:
            raise ExtrapolationError('cell with no division time defined')
        rate, intercept = np.polyfit(times[-npts:], op_values[-npts:], 1)
        val = rate * cell.division_time + intercept

    return y_inv_operator(val)

# %% List of operator acting on timeseries

def _cycle_linear(timeseries):
    times, values = map(np.array, zip(*timeseries))
    a, b = np.polyfit(times, values, 1)
    return a, b


def _cycle_log(timeseries):
    logtimeseries = logarithm(timeseries)
    return _cycle_linear(logtimeseries)


def logarithm(timeseries):
    t, x = map(np.array, zip(*timeseries))
    return (t, np.log(x))


def derivative(timeseries):
    times, values = map(np.array, zip(*timeseries))
    values = np.array(values, dtype='f8')
    delta_t = additive_increments(times)
    delta_v = additive_increments(values)
    newtimes = (times[1:] + times[:-1]) / 2
    der = delta_v / delta_t
    return (newtimes, der)


def logderivative(timeseries):
    times, values = map(np.array, zip(*timeseries))
    values = np.array(values, dtype='f8')
    delta_t = additive_increments(times)
    delta_v = multiplicative_increments(values)
    newtimes = (times[1:] + times[:-1])/2.
    logder = np.log(delta_v) / delta_t
    return (newtimes, logder)


# %% list of operators acting on 1-D arrays

def additive_increments(ar):
    """Computes step-wise additive increments.

    Parameter
    ---------
    ar : 1d Numpy ndarray, n items

    Returns
    -------
    1d Numpy ndarray, n-1 items
    """
    a = ar[1:]
    b = ar[:len(ar)-1]
    return a-b


def multiplicative_increments(ar):
    """Computes step-wise multiplicative increments.

    Parameter
    ---------
    ar : 1d Numpy ndarray, n items

    Returns
    -------
    1d Numpy ndarray, n-1 items
    """
    a = ar[1:]
    b = ar[:len(ar)-1]
    return a/b


# %% functions acting on structured arrays

def compute_secondary_observables(data):
    """Computes secondary observables and extends matrix of observables.

    Argument
    --------
    data -- structured array
        must contains following fields: length, width, fluo, area, time

    Returns
    -------
    out -- structured array
        new fields are added (check `out.dtype.names`)
    """
    ell, w, fluo, area, time = map(np.array,
                                   zip(*data[['length',
                                              'width',
                                              'fluo',
                                              'area',
                                              'time']])
                                   )
    if len(time) > 1:
        delta_t = time[1]-time[0]
        age = (time - time[0] + delta_t/2.)/(time[-1] - time[0] + delta_t)
    else:
        age = np.nan
    volume = spherocylinder_volume(ell, w)
    concentration = fluo/volume
    density = fluo/area
    ALratio = area/ell
    out = append_fields(data,
                        ['volume',
                         'concentration',
                         'density',
                         'ALratio',
                         'age'],
                        [volume,
                         concentration,
                         density,
                         ALratio,
                         age],
                        usemask=False, fill_value=np.nan)
    return out





# %% specific functions

def spherocylinder_volume(length, width):
    """Returns volume of sphero-cylinder.

    Arguments
    ---------
    length -- float
    width -- float

    Returns
    -------
    volume -- float

    Notes
    -----
    sphero-cylinder is composed of a cylinder of height `h`, radius `R`,
    and two hemispheres of radius `R` at each side.

    length is `h + 2*R`
    width is `2*R`
    """
    return np.pi/4.*width**2*(length-width/3.)


def gaussian_smooth(xdata, ydata, sigma=1., x=None):
    """Returns Gaussian smoothed signal.

    Arguments
    ---------
    xdata -- ndarray, co-ordinate of data points
    ydata -- ndarray, ordinates of data points
    sigma -- float (optional, default 1.), x-scale for smoothing
    x (optional) -- float or ndarray,
       values at which the Gaussian smoothed signal is computed
       default xdata

    Returns
    -------
    sequence of (x, gy) co-ordinates of the Gaussian smoothed signal.
    """
    # convert axis for data
    xx = np.expand_dims(xdata, axis=1)
    yy = np.expand_dims(ydata, axis=1)

    def g(t, t0, s):
        return np.exp(-(t-t0)**2/(2.*s**2))/(np.sqrt(2.*np.pi)*s)

    def num(t):
        return np.sum(g(t, xx, sigma)*yy, axis=0)

    def den(t):
        return np.sum(g(t, xx, sigma), axis=0)

    if x is not None:
        u = np.array(x)
    else:
        u = xdata
    return zip(u, num(u)/den(u))


def smooth_timeseries(t, x, option='mean', sigma=1.):
    """Smooth timeseries.

    Arguments
    ---------
    t -- array of timepoints
    x -- array of values to be smoothen, same length as t

    Parameters
    ----------
    option -- str, default 'mean'
        method to choose among ['mean', 'gausssian']
    sigma -- float, default 1
        range of t over which Gaussian smoothing is performed (parameter)

    Returns
    -------
    t, y
        t -- array of times
        y -- array of smoothen values
    """
    nbr = len(t)
    y = np.zeros_like(x)
    if option == 'gaussian':
        return zip(*gaussian_smooth(t, x, sigma=sigma))
    if option == 'mean':
        for k in range(1, nbr-1):
            y[k] = np.mean(x[k-1:k+2])
        y[0] = np.mean(x[:2])
        y[-1] = np.mean(x[-2:])
        return t, y


def get_smooth_cell_timeseries(cell, observable='width'):
    """Returns local 3-points averaged data."""
    time, obs = zip(*cell.data[['time', observable]])
    smobs = []
    if len(time) > 2:
        # average over 3 points within segment
        smobs = [np.mean(obs[i-1:i+2]) for i in range(1, len(obs)-1)]
        # parent
        parent = cell.parent
        if parent and parent.data is not None:
            pt, pobs = zip(*parent.data[['time', observable]])
            smobs.insert(0, np.mean([pobs[-1], obs[0], obs[1]]))
        else:
            # average over first two points
            smobs.insert(0, np.mean(obs[:2]))
        nextobs = []
        for ch in cell.childs:
            chtimes, chobs = zip(*ch.data[['time', observable]])
            nextobs.append(chobs[0])
        if nextobs:
            smobs.append(np.mean([obs[-2], obs[-1], np.mean(nextobs)]))
        else:
            smobs.append(np.mean(obs[-2:]))
    if smobs:
        return zip(time, smobs)
    else:
        return zip(time, obs)


def show_jumps(t, x, threshold=3., mode='multiplicative'):
    """Show data jump over given threshold

    Arguments
    ---------
    t -- ndarray
        measurement times (in hours)
    x -- ndarray (same size as t)
        measurement values

    Parameters
    ----------
    threshold_jump : float
        threshold over which jump is found
        as multiplicative increase x(t+1)/x(t)

    Returns
    -------
    Boolean ndarray with True value for points where there is a jump

    Notes
    -----
    1) These jumps are experimental errors:
       need to be discarded for further analysis.

    2) threshold_jump=1.20 corresponds to 3 doublings/hr
       with one point every 5 mins (2**(3*1/12)=1.19...)
    """
    tdata = np.array(t)
    xdata = np.array(x)
    valids = np.array([False for i in range(len(t))])
    for k in range(len(t)-1):
        try:
            delta_t = tdata[k+1] - tdata[k]
            dot_log = (xdata[k+1]-xdata[k])/(xdata[k] * delta_t)
            if dot_log > threshold*np.log(2.):
                valids[k+1] = True
        except ZeroDivisionError:
            continue
    return valids


if __name__ == '__main__':
    x = np.arange(50, dtype=float)
    y = np.array(len(x) * [np.nan, ])
    anterior_x = np.arange(-20, 0, dtype=int)
    anterior_y = np.array(len(anterior_x) * [np.nan, ])
    anterior_y[np.arange(-20, 0, 5, dtype=int)] = 4.
    y[np.arange(0, len(x), 5, dtype=int)] = 2.
    r, f, xx, yy = compute_rates(x, y, x_break=-.5,
                                 anterior_x=anterior_x, anterior_y=anterior_y,
                                 dt=1, time_window=15.)