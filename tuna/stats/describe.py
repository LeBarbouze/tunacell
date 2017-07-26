#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
tuna package
============

stats module
~~~~~~~~~~~~

describe.py
-----------
"""

from copy import deepcopy

from ..core import Container

from .. import filters
from ..io import fancy


def cell_count(container):
    """Count the number of cells in container.

    Argument
    --------
    container -- Container instance

    Returns
    -------
    counter -- dictionary
       key: number of frames
       value: number od cells having `key` frames
    """
    cells = [node for tree in container.trees for node in tree.all_nodes()]
    counter = {}
    for ce in cells:
        if ce.data is None:
            size = 0
        else:
            size = len(ce.data)
        if size in counter.keys():
            counter[size] += 1
        else:
            counter[size] = 1
    return counter


def tree_count(container):
    """Count the number of trees in container.

    Argument
    --------
    container -- Container instance

    Returns
    -------
    counter -- dictionary
       key: tree depth
       value: number of trees having depth `key`
    """
    counter = {}
    for tr in container.trees:
        size = tr.depth()
        if size in counter.keys():
            counter[size] += 1
        else:
            counter[size] = 1
    return counter


def lineage_count(container):
    """Count the number of independent lineages in container.

    Argument
    --------
    container -- Container instance

    Returns
    -------
    counter -- dictionary
       key: number of cells within lineage
       value: number of lineages having `key` cellsl
    """
    counter = {}
    for tree in container.trees:
        idseqs = tree.decompose()
        for idseq in idseqs:
            size = len(idseq)
            if size in counter.keys():
                counter[size] += 1
            else:
                counter[size] = 1
    return counter


def icount(container, what='cell'):
    if what == 'cell':
        func = cell_count
    elif what == 'tree':
        func = tree_count
    elif what == 'lineage':
        func = lineage_count
    return func(container)


def update_counter(main, add):
    for key, val in add.items():
        if key in main.keys():
            main[key] += val
        else:
            main[key] = val
    return


def count_container(container, filts):
    counter = {'cell': [{} for fi in filts],
               'tree': [{} for fi in filts],
               'lineage': [{} for fi in filts]
               }
    for index, filt in enumerate(filts):
        fcontainer = deepcopy(container)
        fcontainer.postfilter(boofunc=filt, exonerate_root=True)
        for obj, counters in counter.items():
            update_counter(counters[index], icount(fcontainer, what=obj))
        del fcontainer
    return counter


def count_experiment(exp, filts, testing=True):
    """Count lineages and samples.
    """
    if testing:
        upper = 20
        print 'Counting features in testing mode (looping over {} containers)'.format(upper)
    else:
        upper = None
    counter = {'cell': [{} for fi in filts],
               'tree': [{} for fi in filts],
               'lineage': [{} for fi in filts]
               }

    for ifn, container in enumerate(exp.iter_container(size=upper)):
        ccount = count_container(container, filts)
        for obj, counters in ccount.items():
            for index, cts in enumerate(counters):
                update_counter(counter[obj][index], cts)
    return counter


def describe(exp, seqfilters=[], testing=True):
    # setting filter set
    nofilt = filters.cells.FilterCellAny()
    filts = [nofilt] + seqfilters
    if len(seqfilters) > 1:
        combfilt = filters.cells.FilterAND(seqfilters)
        filts += [combfilt]

    counter = count_experiment(exp, filts=filts, testing=testing)

    fancy.pst('Counting cells, trees, lineages')
    for what, count in counter.items():
        store_tot = []
        # printing result
        fancy.title('Counting {}s for each filter'.format(what))
        for index, filt in enumerate(filts):
            with fancy.indent(4):
                fancy.title('{}. filter: {}'.format(index, filt), level=1)
                mat = [['size', 'count']]
                tot = 0
                for key, val in count[index].items():
                    tot += val
                    mat.append('{}\t{}'.format(key, val).split('\t'))
                store_tot.append(tot)
                mat += [[' ', ' '], ['Total', '{}'.format(tot)]]
                fancy.table(mat, header=False, col_sep=1, transpose=True)

        # and the totals table
        with fancy.indent(2):
            fancy.title('total counts of {}s per filter'.format(what), level=2)
            mat = [['filter', 'total count']]
            for index, filt in enumerate(filts):
                mat.append('{}. {}\t{}'.format(index,
                                               repr(filt)[:25],
                                               store_tot[index]
                                               ).split('\t')
                           )
            fancy.table(mat, header=True)
    return
