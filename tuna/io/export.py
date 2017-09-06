#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
tuna package
============

io/export.py module
-------------------
"""
import os
import datetime
from StringIO import StringIO

from tuna.base import datatools
from .. import filters


def output_by_lineages(exp,
                       filtcells=filters.cells.FilterCellAny(),
                       export_path='.',
                       outname=None,
                       testing=True):
    """Output data within lineages, with Marco's requested format.

    Parameters
    -----------
    exp : Experiment instance

    Parameters
    ----------
    filtcells -- FilterCell instance
    filtrees -- FilterTree instance
    filtlin -- FilterLineage instance
    outname -- string, filename for export
    testing -- bool

    TODO: edit README file (automatic), using filter labels
    """
    export_path = os.path.abspath(os.path.expanduser(export_path))

    observables = ['time', 'length', 'width', 'area', 'fluo', 'xCM', 'yCM',
                   'ALratio', 'volume', 'concentration', 'density', 'age']

    count_lineage = 0

    today = datetime.date.today()

    if outname is None:
        bn = 'flat-lineages_experiment-' + exp.metadata.label
    else:
        bn = outname.split('.txt')[0]

    # edit README from 'edit_readme.md'
    readmefilename = os.path.join(export_path, 'readme_' + bn + '.md')
    freadme = open(readmefilename, 'w')

    local_dir = os.path.dirname(__file__)
    print local_dir
    editreadmefilename = os.path.join(local_dir, 'edit_readme.md')
    with open(editreadmefilename, 'r') as editread:
        for line in editread.readlines():
            if line[0] != '!':
                freadme.write(line)
    freadme.close()

    # write proper output
    fname = os.path.join(export_path, bn + '.txt')

    fout = open(fname, 'w')
    fout.write('! HEADER\n')
    metafeed = StringIO(repr(exp.metadata))
    for line in metafeed.readlines():
        fout.write('! ' + line)
    metafeed.close()
    fout.write('! \n')
    fout.write('! Independent lineages -- generated: {}\n'.format(today))
    fout.write('! Filter set:\n')
    filterfeed = StringIO(repr(filtcells))
    for line in filterfeed.readlines():
        fout.write('! \t' + line)
    filterfeed.close()
    fout.write('\n')
    fout.write('! END HEADER\n\n')

    for ifn, container in enumerate(exp.iter_container()):

        if testing and ifn > 5:
            break

        container.postfilter(filt=filtcells)
        fout.write('CONTAINER: {}\n\n'.format(container.label))

        for itree, tree in enumerate(container.trees):

            rootid = tree.root
            idseqs = tree.decompose()
            fout.write('TREE (rootid: {})\n'.format(rootid))

            for index, idseq in enumerate(idseqs):

                count_lineage += 1
                lineID = ','.join(idseq)
                fout.write('LINEAGE: {}\n'.format(lineID))

                data = False
                for icid, cid in enumerate(idseq):

                    # look for Cell instance
                    c = tree.get_node(cid)

                    # test if data (roots without data)
                    if c.data is None:
                        continue

                    data = True
                    alldat = datatools.make_obs(c.data)

                    times, volumes, fluos = zip(*alldat[['time', 'volume', 'fluo']])

                    tau, rate = None, None
                    voli, volf = None, None
                    fi, ff = None, None
                    if c.parent is not None and c.childs and len(times) > 1:
                        tau = times[-1] - times[0] + (times[1]-times[0])
                        rate, intercept = datatools._cycle_log(times, volumes)
                        voli, volf = datatools.extrapolate_endpoints(times,
                                                                     volumes)
                        fi, ff = datatools.extrapolate_endpoints(times, fluos)

                    buff = 'CELL:{}'.format(cid)
                    if tau is not None:
                        buff += ', tau(mins):{}'.format(tau)
                    if rate is not None:
                        buff += ', mu(db/hr):{:.3e}'.format(rate * 60)
                    if voli is not None:
                        buff += ', volume_i(um^3):{:.3e}'.format(voli)
                    if volf is not None:
                        buff += ', volume_f(um^3):{:.3e}'.format(volf)
                    if fi is not None:
                        buff += ', fluo_i(a.u.):{:.3e}'.format(fi)
                    if ff is not None:
                        buff += ', fluo_f(a.u.):{:.3e}'.format(ff)
                    buff += '\n'

                    fout.write(buff)

                    for obs in observables:
                        values = alldat[obs]
                        valstr = obs + '\t' + '\t'.join(
                               ['{:.3e}'.format(val) for val in values]) + '\n'
                        fout.write(valstr)
                if data:
                    fout.write('\n')
    fout.close()
    return