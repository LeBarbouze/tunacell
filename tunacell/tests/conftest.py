#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Test configurations, include fixtures"""
import os
import pytest
import shutil

import tunacell
from tunacell.base.experiment import Experiment

path_data = os.path.join(os.path.dirname(tunacell.__file__), 'data')
path_fake_exp = os.path.join(path_data, 'fake')


@pytest.fixture(scope='session')
def fake_exp():
    exp = Experiment(path_fake_exp, count_items=True)
    yield exp
    # remove analysis path once tests have been performed
    analysis_path = os.path.join(path_fake_exp, 'analysis')
    if os.path.exists(analysis_path):
        shutil.rmtree(analysis_path)
    # remove internals path
    if exp.path_internals.exists():
        shutil.rmtree(str(exp.path_internals))