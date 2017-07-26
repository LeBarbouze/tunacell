#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
tuna package
============

io module
~~~~~~~~~
Input/Output stuff

fancy.py
--------
where we use clint package to fancy standard output
"""

from clint.textui import puts, indent, colored
from clint.textui import columns

import StringIO


class FancyPrintingError(Exception):
    pass


class FancyTableError(FancyPrintingError):
    pass


def table(mat, header=False, col_sep=3, indent_table=1, transpose=False):
    """Print fancy table

    Argument
    --------
    mat -- list of sequences of strings

        >>> mat = [['animal', 'couleur', 'age'],
                   ['chien', 'sable', '3'],
                   ['chat', 'noir', '7']]
    """
    mat_shape = len(mat), len(mat[0])
    if transpose:
        mat = [[mat[i][j] for i in range(mat_shape[0])]
                for j in range(mat_shape[1])]
    nbr_cols = len(mat[0])
    for line in mat:
        if len(line) != nbr_cols:
            raise FancyTableError('Unadequate number of columns')
    # detect maximum size for each column
    size_cols = nbr_cols * [0]
    for line in mat:
        for index, item in enumerate(line):
            if len(item) > size_cols[index]:
                size_cols[index] = len(item)
    # right justify elements, line per line
    wmat = []
    for index, line in enumerate(mat):
        justified_line = []
        for j, item in enumerate(line):
            justified_line.append(item.rjust(size_cols[j]))
        wmat.append(justified_line)
    # reajust col size given col_sep
    size_cols = [sc + col_sep - 1 for sc in size_cols]
    # write table
    with indent(indent_table, quote=''):
        start = 0
        if header:
            elems = tuple(map(list, zip(wmat[0], size_cols)))
            puts(columns(*elems))
            start = 1
            underlines = tuple(map(list,
                                   zip([len(item)*'-' for item in wmat[0]],
                                        size_cols)
                                   )
                              )
            puts(columns(*underlines))
        for line in wmat[start:]:
            elems = tuple(map(list, zip(line, size_cols)))
            puts(columns(*elems))
        puts()  # blank line at the end
    return


def title(string, level=0, indent_newlines=0):
    items = string.split('\n')
    size = max(map(len, items))
    # set marker
    if level == 0:
        marker = '='
    elif level == 1:
        marker = '*'
    elif level == 2:
        marker = '~'
    else:
        marker = '-'
    # printing first line
    puts(items[0])
    # test whether new lines to print
    if len(items) > 1:
        with indent(indent_newlines + 1, quote=''):
            for newline in items[1:]:
                puts(newline)
        size = size + indent_newlines + 1
    puts(size * marker)
    puts()
    return


def pst(string):
    "Print Super Title"
    items = string.split('\n')
    size = max(map(len, items))
    marker = '='
    puts()  # one blank line before
    with indent(12, quote=''):

        puts((size + 8) * marker)
        for line in items:
            puts('|   ' + line + '   |')
        puts((size + 8) * marker)
    puts()  # one blank line after
    return


def tree(tree, indent_tree=2):
    "Print Trees (tree must be treelib.Tree instance"
    ft = StringIO.StringIO(str(tree))
    with indent(indent_tree, quote=''):
        for line in ft.readlines():
            puts(line[:-1])
    puts()
    return