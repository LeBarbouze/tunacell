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
            msg = ('Arrays of different lengths! Check x and y input\n'
                   'x : {}'.format(x) + '\n'
                   'y : {}'.format(y))
            raise ValueError(msg)
        self._x = np.array(x)
        self._y = np.array(y)
        self.valid = np.where(np.logical_not(np.isnan(y)))
        return

    @property
    def x(self):
        return self._x

    @property
    def y(self):
        return self._y

    @property
    def clear_x(self):
        """Returns x coordinates for which y is valid"""
        if len(self._y) > 0:
            return self._x[self.valid]
        else:
            return np.array([], dtype=float)

    @property
    def clear_y(self):
        """Returns y coordinates when y is valid"""
        if len(self._y) > 0:
            return self._y[self.valid]
        else:
            return np.array([], dtype=float)

    @property
    def clear(self):
        """Returns coordinates cleared off NaNs"""
        return Coordinates(self.clear_x, self.clear_y)

    def as_array(self, x_name=None, y_name=None):
        if x_name is None:
            x_name = 'x'
        if y_name is None:
            y_name = 'y'
        array = np.array(zip(self.clear_x, self.clear_y),
                         dtype=[(x_name, 'f8'), (y_name, 'f8')])
        return array


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
    anterior_rates : 1d ndarray
        array of rates estimated at anterior_x array
    anterior_fits : 1d ndarray
        array of fitter values at anterior_x array
    used_x : 1d ndarray
        concatenation with derivative continuity hypothesis of anterior_x and x
    used_y : 1d ndarray
        concatenation with derivative continuity hypothesis of anterior_y and y
    """
    # check array lengths and associate cleaning ordinate NaNs
    coords = Coordinates(x, y)
    nans_coords = len(x) * [np.nan, ]
    anteriors = Coordinates(anterior_x, anterior_y)
    nans_anteriors = len(anterior_x) * [np.nan, ]

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

    # find period (can be different from dt when multiple acquisition periods)
    no_data = len(coords.clear_x) == 0
    too_few_data = coords.clear_x[-1] - coords.clear_x[0] < time_window
    if no_data or too_few_data:
        return nans_coords, nans_coords, nans_anteriors, nans_anteriors, [], []

    # coords.clear_x is necesarily of length >=2
    period = np.amin(coords.clear_x[1:] - coords.clear_x[:-1])

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

    op_y = y_operator(coords.clear_y)

    if len(coords.clear_x) >= join_points:
        # fit to at least join_points, more if possible
        r, i = np.polyfit(coords.clear_x[:n_joints], op_y[:n_joints], 1)
        op_y_break = i + r * x_break

    if testing:
        msg = ('Data to fit:\n'
               'x : {}'.format(coords.clear_x) + '\n'
               'y : {}'.format(coords.clear_y))
        print(msg)

    # try to use anterior values : compute break offset
    trans_op_ay = []  # translated, operated anterior values; default: empty
    offset = None
    if op_y_break is not None and len(anteriors.clear_x) > 0:
        op_ay = y_operator(anteriors.clear_y)
        # 2 checks:
        #   1. there at enough points to get the final value estimate
        cdt1 = len(anteriors.clear_x) >= join_points
        #   2. initial value is determined
        cdt2 = op_y_break is not None
        if cdt1 and cdt2:
            r, i = np.polyfit(anteriors.clear_x[-n_joints:],
                              op_ay[-n_joints:], 1)
            op_ay_break = i + r * x_break
            if testing:
                msg = ('Extrapolated anterior value at break:\n '
                       'break time, value: '
                       '{}, {}'.format(x_break, y_inv_operator(op_ay_break)))
                print(msg)
                msg = ('Anterior data used for fitting:\n'
                       'x : {}'.format(anteriors.clear_x) + '\n'
                       'y : {}'.format(anteriors.clear_y))
                print(msg)
            # adjust values by translating
            offset = op_ay_break - op_y_break
            trans_op_ay = op_ay - offset
    if testing:
        msg = ('Translated, operated anterior values are:\n'
               'y\' {}'.format(trans_op_ay))
        print(msg)
        print()

    # concatenate operated values if enough previous values to comppute offset
    if len(trans_op_ay) > 0:
        all_x = np.concatenate([anteriors.clear_x, coords.clear_x])
        all_op_y = np.concatenate([trans_op_ay, op_y])
    else:
        all_x = coords.clear_x
        all_op_y = op_y
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
            if testing:
                msg = ('+ window range does not intercept x range: next')
                print(msg)
            continue
        # reduce to points to fit
        lower = all_x > t_start
        upper = all_x <= t_stop
        boo = np.logical_and(lower, upper)
        if testing:
            print('+ window:', end=' ')
            print('{} < t < {} '.format(t_start, t_stop), end=' ')
            print('({} points)'.format(len(all_x[boo])))
        # check if n_points > 3
        if len(all_x[boo]) < n_points:
            if testing:
                print('Not enough points on this time window')
                print('({} instead of {})'.format(len(all_x[boo]), n_points))
                print('next')
            # not enough point in this window: insert NaN
            rate_op_y[index] = np.nan
            fit_op_y[index] = np.nan
        else:
            rate, intercept = np.polyfit(all_x[boo], all_op_y[boo], 1)
            fit_op_y[index] = rate * time_eval + intercept
            rate_op_y[index] = rate
            if testing:
                msg = ('time          : {}'.format(time_eval) + '\n'
                       'fitted value  : {}'.format(fit_op_y[index]) + '\n'
                       'computed rate : {}'.format(rate))
                print(msg)

    # new values: linear interpolation
    fit_op_coords = Coordinates(fit_x, fit_op_y)
    f = interp1d(fit_op_coords.clear_x, fit_op_coords.clear_y, kind='linear',
                 assume_sorted=True, bounds_error=False)
    # interpolation may be defined over both x range and anterior_x range
    out_y = np.array(len(coords.x) * [np.nan, ])  # initialize to NaNs
    out_anterior_y = np.array(len(anteriors.x) * [np.nan, ])
    # valid data points: interpolation at initial x coordinates
    out_y[coords.valid] = y_inv_operator(f(coords.clear_x))
    if len(anteriors.valid) > 0 and offset is not None:
        out_anterior_y[anteriors.valid] = y_inv_operator(f(anteriors.clear_x) + offset)

    # rates : linear interpolation
    rate_coords = Coordinates(fit_x, rate_op_y)
    f = interp1d(rate_coords.clear_x, rate_coords.clear_y, kind='linear',
                 assume_sorted=True, bounds_error=False)
    out_rate = np.array(len(coords.x) * [np.nan, ])
    out_anterior_rate = np.array(len(anteriors.x) * [np.nan, ])
    # valid data points: interpolation at initial x coordinates
    out_rate[coords.valid] = f(coords.clear_x)
    if len(anteriors.valid) > 0 and offset is not None:
        out_anterior_rate[anteriors.valid] = f(anteriors.clear_x)

    return out_rate, out_y, out_anterior_rate, out_anterior_y, all_x, all_y


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

# %% List of operator acting on Coordinates


def _cycle_linear(coords):
    if len(coords.valid) < 2:
        raise ValueError('Not enough valid coordinates')
    a, b = np.polyfit(coords.clear_x, coords.clear_y, 1)
    return a, b


def _cycle_log(coords):
    logcoords = logarithm(coords)
    return _cycle_linear(logcoords)


def logarithm(coords):
    if len(coords.valid) < 2:
        raise ValueError('Not enough valid coordinates')
    return Coordinates(coords.x, np.log(coords.y))


def derivative(coords):
    delta_x = additive_increments(coords.clear_x)
    delta_y = additive_increments(coords.clear_y)
    new_x = (coords.clear_x[1:] + coords.clear_y[:-1])/2.
    new_y = delta_y/delta_x
    # interpolate to associate to initial times
    f = interp1d(new_x, new_y, kind='linear', assume_sorted=True, bounds_error=False)
    out_y = np.array(len(coords.x) * [np.nan, ])
    if len(coords.valid) > 0:
        out_y[coords.valid] = f(coords.clear_x)
    return Coordinates(coords.x, out_y)


def logderivative(coords):
    delta_x = additive_increments(coords.clear_x)
    delta_y = multiplicative_increments(coords.clear_y)
    new_x = (coords.clear_x[1:] + coords.clear_y[:-1])/2.
    new_y = np.log(delta_y)/delta_x
    # interpolate to associate to initial times
    f = interp1d(new_x, new_y, kind='linear', assume_sorted=True, bounds_error=False)
    out_y = np.array(len(coords.x) * [np.nan, ])
    if len(coords.valid) > 0:
        out_y[coords.valid] = f(coords.clear_x)
    return Coordinates(coords.x, out_y)


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
    r, f, ar, af, xx, yy = compute_rates(x, y, x_break=-.5, scale='linear',
                                         anterior_x=anterior_x,
                                         anterior_y=anterior_y,
                                         dt=1, time_window=15., testing=True)
    coords = Coordinates(np.concatenate([anterior_x, x]),
                         np.concatenate([af, f]))
    array = coords.as_array(x_name='time', y_name='value')
    print(array['time'])
    print(array['value'])
    