#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

"""
from __future__ import print_function

import pytest
import numpy as np
import itertools

from tunacell.simu.base import SimuParams, DivisionParams, SampleInitialSize
from tunacell.simu.ou import OUSimulation


@pytest.fixture
def simu_params():
    sp = SimuParams(nbr_container=5, nbr_colony_per_container=1,
                    start=0., stop=400., period=5., seed=42)
    return sp


def test_division_params():
    """Test DivisionParams"""
    lambdas = np.linspace(0, 1, 5)  # 5 values
    sizes = np.linspace(1, 5, 5)
    modes = ['none', 'gamma']
    flucts = np.linspace(0.1, 1, 5)

    growth_rates = np.log(2.) / (60. * np.array([20., 33.3, 47.7, 60., 1001.]))

    for item in itertools.product(lambdas, sizes, modes, flucts):
        div_lambda, div_size, div_mode, div_sd_to_mean = item
        dp = DivisionParams(div_lambda, div_size, div_mode, div_sd_to_mean)

        # avoid problem by reaching division size cutoff
        birth_sizes = np.linspace(div_size/10., div_size, 5, endpoint=False)

        for bs, mu in itertools.product(birth_sizes, growth_rates):
            assert dp.rv(birth_size=bs, growth_rate=mu) > 0.


def test_sample_initial_size():
    """Test SampleInitialSize"""
    sizes = np.linspace(1, 10, 7)
    # mode: fixed
    for size in sizes:
        bs = SampleInitialSize(birth_size_mean=size)
        assert bs.rv() == size
    # mode: lognormal
    flucts = np.linspace(0.1, 3, 7)
    for size, fluct in itertools.product(sizes, flucts):
        bs = SampleInitialSize(birth_size_mean=size,
                               birth_size_mode='lognormal',
                               birth_size_sd_to_mean=fluct)
        assert bs.rv() > 0.


def test_ou_simulation_container_iteration():
    ou = OUSimulation()
    expected_nbr_container = ou.simuParams.nbr_container
    counter = 0
    for container in ou.iter_containers():
        counter += 1
        assert container.filetype == "simulations"
    assert counter == expected_nbr_container


