#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Testing container features
"""
from __future__ import print_function

import pytest
import tuna
import os

from tuna.base.experiment import Experiment

path_data = os.path.join(os.path.dirname(tuna.__file__), 'data')
path_fake_exp = os.path.join(path_data, 'fake')


@pytest.fixture(scope="module")
def container():
    """Get first container"""
    exp = Experiment(path_fake_exp)
    return exp.get_container('container_01')


def test_container_data(container):
    assert container.data['cellID'][0] ==  1
    assert container.data['parentID'][0] == 0
    assert container.data['time'][0] == 0.
    assert container.data['value'][0] == 1.
    

