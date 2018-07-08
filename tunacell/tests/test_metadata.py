#!/usr/bin/env python
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
import yaml

import tunacell
#from tunacell.base.experiment import Experiment
from tunacell.io.metadata import Metadata, LocalMetadata, load_from_yaml, load_from_csv

path_data = os.path.join(os.path.dirname(tunacell.__file__), 'data')
path_fake_exp = os.path.join(path_data, 'fake')

path_yaml_examples = os.path.dirname(os.path.abspath(__file__))

#@pytest.fixture(scope='module')
#def exp():
#    expe = Experiment(path_fake_exp)
#    return expe
#
#
#def test_experiment_meta(exp):
#    assert exp.label in exp.metadata.index
#    row = exp.metadata.loc[exp.label]
#    assert row.author == 'Joachim Rambeau'
#    assert row.period == 3.14
#    assert row.strain == 'Thunnus alalunga'
#
#
#def test_container_meta(exp):
#    assert 'container_02' in exp.metadata.index
#    row = exp.metadata.loc['container_02']
#    assert row.author == 'Joachim Rambeau'
#    assert row.period == 3.14
#    assert row.strain == 'Thunnus thynnus'


@pytest.fixture
def simple():
    path = os.path.join(path_yaml_examples, 'simple_metadata.yml')
    stream = open(path, 'r')
    yield yaml.load_all(stream)
    stream.close()


@pytest.fixture
def layers():
    path = os.path.join(path_yaml_examples, 'layers_metadata.yml')
    stream = open(path, 'r')
    yield yaml.load_all(stream)
    stream.close()


@pytest.fixture
def from_yaml():
    return load_from_yaml(os.path.join(path_fake_exp, 'metadata.yml'))


@pytest.fixture
def from_csv():
    return load_from_csv(os.path.join(path_fake_exp, 'metadata.csv'), sep=',')


def test_simple_load(simple):
    """Checks that constructor works and local metadata is assigned to top"""
    md = Metadata(simple)
    assert isinstance(md, Metadata)
    assert isinstance(md.experiment, LocalMetadata)


def test_simple_content(simple):
    md = Metadata(simple)
    assert md['experiment'] == 'fake'
    assert md['period'] == 5.0
    assert md['author'] == 'J. Rambeau'
    assert md['species'] == 'e.coli'
    assert md['temperature'] == {'value': 37, 'units': 'C'}


def test_simple_period(simple):
    md = Metadata(simple)
    assert md.period == 5.0
    assert md.experiment.period == 5.0


def test_layers_load(layers):
    md = Metadata(layers)
    assert len(md.locals) == 3  # 3 containers


def test_layers_period(layers):
    md = Metadata(layers)
    assert md.period == 5.0  # should pick minimum


def test_layers_simple_container(layers):
    md = Metadata(layers)
    loc = md.from_container('container_007')
    assert loc['species'] == 'james bond'
    assert md['species'] == 'e.coli'


def test_layers_nested_containers(layers):
    md = Metadata(layers)
    loc = md.from_container('container_xy')
    assert loc['medium'] == 'M9Fru'
    assert loc['species'] == 'jellyfish'
    assert loc['period-ch2'] == 20.0
