#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
test_cell
"""
from __future__ import print_function

import pytest
import numpy as np

from tunacell.base.cell import Cell, filiate_from_bpointer
from tunacell.base.colony import Colony, build_recursively_from_cells


## Fixture for cell

@pytest.fixture
def cells():
    """Builds a list of 10 cells with ids from '0' to '9'"""
    cells = [Cell(identifier='0'), ]
    for index in range(1, 10):
        cells.append(Cell(identifier='{}'.format(index)))
    return cells


@pytest.fixture
def binary_division_cells(cells):
    # we want to label cells sequentially given their depth in filiation
    # e.g. [[0], [1 2], [3 4 5 6], ...]
    # first define how to slice an array of [0 1 2 3 4 5 6 ...]
    sls = []
    for k in range(1, 5):
        sls.append(slice(sum((2**i for i in range(k-1))), sum((2**i for i in range(k)))))
    # second, slice the array of 10 cells
    ll = [cells[sl] for sl in sls]

    # third associate parent
    for irow, row in enumerate(ll[1:], start=1):
        for icol, col in enumerate(row):
            col.bpointer = ll[irow - 1][icol // 2].identifier
    return cells


## Fixtures for colonies

@pytest.fixture
def tree():
    """Example from `treelib`_
    
    .. _treelib: https://github.com/caesar0301/treelib/blob/master/examples/family_tree.py
    """
    tree = Colony()
    tree.create_node("Harry", "harry")  # root node
    tree.create_node("Jane", "jane", parent="harry")
    tree.create_node("Bill", "bill", parent="harry")
    tree.create_node("Diane", "diane", parent="jane")
    tree.create_node("Mary", "mary", parent="diane")
    tree.create_node("Mark", "mark", parent="jane")
    return tree

## Test functions for cell

def test_cell_parent(cells):
    root = cells[0]
    for index in range(1, 10):
        cell = cells[index]
        cell.parent = root
    for index in range(1, 10):
        assert cell.parent == root
    

def test_cell_childs(cells):
    root = cells[0]
    for index in range(1, 10):
        cell = cells[index]
        root.childs = cell
    assert root.childs == cells[1:]
    

def test_cell_division_timing(cells):
    root = cells[0]
    couples = [(int(root.identifier), t) for t in [0, 1, 2]]
    root.data = np.array(couples, dtype=[('cellID', 'u2',), ('time', 'f8')])
    for index in range(1, 10):
        cell = cells[index]
        times = [4, 5, 6]
        couples = [(cell.identifier, t) for t in times]
        cell.data = np.array(couples, dtype=[('cellID', 'u2',), ('time', 'f8')])
        cell.parent = root
        cell.set_division_event()
    
    assert root.division_time == 3
    
    for index in range(1, 10):
        assert cell.birth_time == 3
    
    
def test_cell_build_filiation(binary_division_cells):
    cells = binary_division_cells
     # make filiation in place
    filiate_from_bpointer(cells)
    
    assert cells[1] in cells[0].childs
    assert cells[2] in cells[0].childs
    assert cells[1].parent == cells[0]
    assert cells[2].parent == cells[0]
    assert cells[3] in cells[1].childs
    assert cells[4] in cells[1].childs
    assert cells[3].parent == cells[1]
    assert cells[4].parent == cells[1]
    assert cells[5] in cells[2].childs
    assert cells[6] in cells[2].childs
    assert cells[5].parent == cells[2]
    assert cells[6].parent == cells[2]
    
    
def test_colony_full_decomposition(tree):
    tree.decompose(independent=False, seed=42)
    assert tree._decomposition['independent'] is False
    assert tree._decomposition['seed'] == 42
    assert tree.idseqs == [['harry', 'bill'],
                           ['harry', 'jane', 'diane', 'mary'],
                           ['harry', 'jane', 'mark']]


def test_colony_independent_decomposition(tree):
    tree.decompose(independent=True, seed=42)
    assert tree._decomposition['independent'] is True
    assert tree._decomposition['seed'] == 42
    assert tree.idseqs == [['harry', 'jane', 'mark'], ['diane', 'mary'], ['bill']]


def test_colony_recursive_constructor(binary_division_cells):
    """This is effectively used in Container class"""
    cells = binary_division_cells
    filiate_from_bpointer(cells)
    colonies = build_recursively_from_cells(cells)
    assert len(colonies) == 1
    colony = colonies[0]
    assert colony.depth() == 3
    assert colony.level('9') == 3

