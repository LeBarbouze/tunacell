#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
tunacell package
============

test suite
~~~~~~~~~~~~~~~~~~~~~~~~~
"""
from __future__ import print_function

import pytest
import os

import tunacell
from tunacell.base.experiment import Experiment


path_data = os.path.join(os.path.dirname(tunacell.__file__), 'data')
path_fake_exp = os.path.join(path_data, 'fake')


@pytest.fixture(scope='module')
def exp():
    expe = Experiment(path_fake_exp)
    return expe


def test_experiment_meta(exp):
    assert exp.label in exp.metadata.index
    row = exp.metadata.loc[exp.label]
    assert row.author == 'Joachim Rambeau'
    assert row.period == 3.14
    assert row.strain == 'Thunnus alalunga'


def test_container_meta(exp):
    assert 'container_02' in exp.metadata.index
    row = exp.metadata.loc['container_02']
    assert row.author == 'Joachim Rambeau'
    assert row.period == 3.14
    assert row.strain == 'Thunnus thynnus'
