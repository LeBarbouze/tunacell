#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Testing filtering features.
"""
from __future__ import print_function

import pytest

from tunacell import Observable

from tunacell.filters.main import FilterSet
from tunacell.filters.cells import (FilterCellAny,
                                FilterCellIDparity,
                                FilterObservableBound
                                )
from tunacell.filters.lineages import FilterLineageWithCellProperty


def test_cell_any():
    filt = FilterCellAny()
    assert filt._type == 'CELL'
    assert filt.label == 'CELL, Always True'
    assert filt._obs == []
    assert filt('whatever')  # always True


def test_unroll():
    o1 = Observable(name='test-obs-1-used-for-cells', raw='length')
    o2 = Observable(name='test-obs-2-used-for-lineages', raw='width')
    f1 = FilterObservableBound(obs=o1, lower_bound=1, upper_bound=2)  # cells
    f2 = FilterObservableBound(obs=o2, lower_bound=3, upper_bound=4)  # cells
    lf = FilterLineageWithCellProperty(cell_filter=f2)  # lineages
    fset = FilterSet(label='new-idea', filtercell=f1, filterlineage=lf)
    assert set(fset.obs) == set([o1, o2])
    
