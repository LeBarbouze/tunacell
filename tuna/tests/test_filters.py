#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Testing filtering features.
"""
from __future__ import print_function

import pytest


from tuna.filters.cells import (FilterCellAny,
                                FilterCellIDparity,
                                )


def test_cell_any():
    filt = FilterCellAny()
    assert filt._type == 'CELL'
    assert filt.label == 'CELL, Always True'
    assert filt._obs == []
    assert filt('whatever')  # always True
