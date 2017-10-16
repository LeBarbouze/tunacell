#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Testing Experiment features.
"""
from __future__ import print_function

import pytest
import os

import tunacell
from tunacell.base.experiment import Experiment
from tunacell.base.container import Container

path_data = os.path.join(os.path.dirname(tunacell.__file__), 'data')
path_fake_exp = os.path.join(path_data, 'fake')


@pytest.fixture(scope='module')
def fake_exp():
    exp = Experiment(path_fake_exp)
    return exp


def test_load_experiment(fake_exp):
    assert isinstance(fake_exp, Experiment)


def test_experiment_attributes(fake_exp):
    attrs = ['label', 'abspath', 'metadata', 'datatype', 'filetype',
             'containers']
    for attr in attrs:
        assert hasattr(fake_exp, attr)


def test_experiment_label(fake_exp):
    assert fake_exp.label == 'fake'


def test_experiment_abspath(fake_exp):
    assert fake_exp.abspath == path_fake_exp


def test_experiment_number_container(fake_exp):
    assert len(fake_exp.containers) == 3


def test_experiment_metadata(fake_exp):
    assert fake_exp.metadata.loc[fake_exp.label, 'author'] == 'Joachim Rambeau'


def test_experiment_datatype(fake_exp):
    assert len(fake_exp.datatype) == 4
    assert 'cellID' in list(zip(*fake_exp.datatype))[0]


def test_experiment_iter(fake_exp):
    meths = ['iter_containers']
    for meth in meths:
        method = getattr(fake_exp, meth)
        assert callable(method)


def test_experiment_get_container(fake_exp):
    container = fake_exp.get_container('container_01')
    assert isinstance(container, Container)
