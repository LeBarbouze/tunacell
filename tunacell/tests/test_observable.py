#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Testing Observable features.
"""
from __future__ import print_function

import pytest
from itertools import product

from tunacell.base.observable import (Observable, ObservableStringError,
                                  FunctionalObservable, set_observable_list,
                                  unroll_func_obs, unroll_raw_obs)


t_params = (('name', ['mysoup', 'donald']),
            ('raw', ['soup', 'duck']),
            ('scale', ['linear', 'log']),
            ('mode', ['dynamics',
                      'birth', 'division',
                      'net-increase-additive', 'net-increase-multiplicative',
                      'rate', 'average']),
            ('join_points', [3, 10, 100]),
            ('timing', ['t', 'b', 'd', 'g']),
            ('tref', [None, 0., 30, 300.]),
            ('differentiate', [False, True]),
            ('local_fit', [False, True]),
            ('time_window', [0., 3.14, 15.])
            )


#@pytest.fixture
#def random_strings():
#    average_string_length = 100
#    number_of_samples = 20
#    lengths = np.random.poisson(lam=average_string_length,
#                                size=number_of_samples)
#
#    def gen_random_strings(sizes):
#        chars = list(string.printable)
#        for size in sizes:
#            array = np.random.randint(0, high=len(chars), size=size)
#            yield ''.join([chars[k] for k in array])
#
#    return gen_random_strings(lengths)


@pytest.fixture(scope='module')
def all_params():
    """Fixture to return an iterator among all parameter sets"""
    def iter_params(params):
        keys, values = zip(*params)
        for items in product(*values):
            yield {key: value for key, value in zip(keys, items)}

    return iter_params(t_params)


def test_observable_init(all_params):
    """Test Observable initialization"""
    for kwargs in all_params:
        obs = Observable(**kwargs)
        for attr in obs._attr_names:
            assert getattr(obs, attr) == kwargs[attr]


def test_observable_str(all_params):
    """This test checks codestring obtained using __str__

    Randomly instantiate a given observable <obs>
    Check that second instance made from str(obs) is identical.
    """
    for kwargs in all_params:
        obs = Observable(**kwargs)
        nobs = Observable(from_string=str(obs))
        for attr in obs._attr_names:
            if attr == 'time_window' and not obs.local_fit:
                # check only if local_fit is True,
                # since when it's false, time_window is not printed in str
                continue
            assert hasattr(nobs, attr)
            assert getattr(nobs, attr) == getattr(obs, attr)


def test_observable_repr(all_params):
    """This test checks string representation obtained using __repr__

    Randomly instantiate a given observable <obs>
    Check that second instance made from str(obs) is identical.
    """
    for kwargs in all_params:
        obs = Observable(**kwargs)
        nobs = eval(repr(obs))
        for attr in obs._attr_names:
            assert hasattr(nobs, attr)
            assert getattr(nobs, attr) == getattr(obs, attr)


def test_unrolling():
    length = Observable(name='length')
    width = Observable(name='width')
    area = FunctionalObservable(name='square-area', f=lambda x, y: x*y, observables=[length, width])
    vol = FunctionalObservable(name='volume', f=lambda x, y: x*y, observables=[area, width])
    raw_obs, func_obs = set_observable_list(vol, filters=[])
    assert length in raw_obs
    assert length not in func_obs
    assert width in raw_obs
    assert width not in func_obs
    assert area in func_obs
    assert area not in raw_obs
    assert vol in func_obs
    assert vol not in raw_obs

def test_unroll_func():
    length = Observable(name='length')
    width = Observable(name='width')
    area = FunctionalObservable(name='square-area', f=lambda x, y: x*y, observables=[length, width])
    vol = FunctionalObservable(name='volume', f=lambda x, y: x*y, observables=[area, width])
    func_obs = list(unroll_func_obs(vol))
    assert func_obs == [area, vol]  # oredering matters for functionalObs
    func_obs = list(unroll_func_obs([vol, area]))
    assert func_obs == [area, vol, area]
