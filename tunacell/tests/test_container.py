#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Testing container features
"""
from __future__ import print_function

import pytest
import tunacell
import os

import sys
if sys.version_info[0] < 3:
    import pathlib2 as pathlib
else:
    import pathlib

from tunacell.base.experiment import Experiment

path_data = os.path.join(os.path.dirname(tunacell.__file__), 'data')
path_fake_exp = os.path.join(path_data, 'fake')


@pytest.fixture(scope="module")
def container():
    """Get first container"""
    exp = Experiment(path_fake_exp)
    return exp.get_container('container_01')


def test_attributes(container):
    if isinstance(container.abspath, pathlib.Path):
        assert str(container.abspath.absolute()) == os.path.join(os.path.abspath(path_fake_exp),
                                              'containers',
                                              'container_01.txt')
    else:
        assert container.abspath == os.path.join(os.path.abspath(path_fake_exp),
                                              'containers',
                                              'container_01.txt')
    assert container.filetype == 'text'
    assert container.label == 'container_01'
    assert container.period == 3.14
    # metadata
    md = container.metadata
    assert md['author'] == 'Joachim Rambeau'
    assert md['strain'] == 'Thunnus alalunga'


def test_container_data(container):
    assert container.data['cellID'][0] == 1
    assert container.data['parentID'][0] == 0
    assert container.data['time'][0] == 0.
    assert container.data['value'][0] == 1.


def test_container_content(container):
    assert len(container.cells) == 6
    assert len(container.trees) == 1
    colony = container.get_colony(2)
    assert colony.root == 1  # root cell
