#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Modules that defines import from/export to HDF5 files

TODO: re-implement
"""
from __future__ import print_function


class H5ParsingError(Exception):
    pass


def get_array(table):
    arr = table.read()
    return arr


#def build_cells(table, container=None, report_NaNs=True):
#    """Read and store Cell instances from table..
#
#    Argument
#    --------
#    table -- tables.Table instance, where data is stored
#
#    Returns
#    -------
#    list of Cell instances
#       Information is stored in:
#           * bpointer: backwards pointer, to parent cell
#           * .data: data as structured array
#    """
#    # 0. get observable names
#    observables = tuple(table.colnames)
#    # 1. find separation in table between cells
#    marks = [0]
#    for row in table.iterrows(start=0, stop=1):
#        cid = row['cellID']  # raw cid, will convert to string below
#    for index, row in enumerate(table.iterrows(start=1), start=1):
#        if row['cellID'] != cid:
#            marks.append(index)
#            cid = row['cellID']
#    sls = [slice(marks[i], marks[i+1])
#           for i in range(len(marks)-1)] + [slice(marks[-1], None)]
#    cells = []
#    errors = []
#    for sl in sls:
#        arr = table.read(start=sl.start, stop=sl.stop)
#        cid, pid = map(str, arr[['cellID', 'parentID']][0])
#        cell = Cell(identifier=cid, container=container)
#        if pid != '0':
#            cell.bpointer = pid
#        cell.data = arr
#        cells.append(cell)
#                # record if NaN values appear
#        if report_NaNs:
#            for label, (dtype, offset) in arr.dtype.fields.items():
#                # NaNs are implemented as np.nan for float types,
#                if 'f' in dtype.kind:
#                    if np.isnan(arr[label]).any():
#                        row = arr[np.isnan(arr[label])]
#                        err = CellNaNError(identifier=cid, row=row,
#                                           container=container, column_label=label,
#                                           column_type=dtype)
#                        errors.append(err)
#                # for integer types, they seem to be replaced by largest value
#                elif ('u' in dtype.kind) or ('i' in dtype.kind):
#                    if np.amax(arr[label]) == np.iinfo(dtype).max:
#                        row = arr[arr[label] == np.iinfo(dtype).max]
#                        err = CellNaNError(identifier=cid, row=row,
#                                           container=container, column_label=label,
#                                           column_type=dtype)
#                        errors.append(err)
#
#    return cells, observables, errors
