#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Testing io.text module.
"""
from __future__ import print_function

import pytest
import os

import tunacell
from tunacell.io import text


path_data = os.path.join(os.path.dirname(tunacell.__file__), 'data')
path_fake_exp = os.path.join(path_data, 'fake')


@pytest.fixture
def datatype():
    return text.datatype_parser(os.path.join(path_fake_exp, 'descriptor.csv'))


@pytest.fixture
def fname():
    return os.path.join(path_fake_exp, 'lineages', 'container_01.txt')


def test_check_up(fname):
    path = os.path.dirname(fname)
    fn = text._check_up('metadata.csv', path, level=1)
    assert os.path.basename(fn) == 'metadata.csv'
    with pytest.raises(text.MissingFileError):
        _ = text._check_up('metadata.csv', path, level=0)


def test_filename_parser(fname):
    containers = text.find_containers(path_fake_exp)
    assert 'container_01.txt' in list(map(lambda item: item.name, containers))
    assert 'container_02.txt' in list(map(lambda item: item.name, containers))
    assert 'container_03.txt' in list(map(lambda item: item.name, containers))
    assert len(containers) == 4
