#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
    supersegger_trials
    ^^^^^^^^^^^^^^^^^^

    Example to test/use supersegger files
"""

from tunacell import Experiment, Parser, Observable
from tunacell.plotting.samples import SamplePlot


def iterate(exp):
    count = 0
    for cell in exp.iter_cells():
        count += 1
    print('{} cells'.format(count))
    print('last: {}'.format(cell))
    return cell

def plot_samples(exp):
    print(isinstance(exp, Experiment))
    parser = Parser(exp)
    parser.add_sample(1)  # random sample
    colony = parser.get_colony(0)
    print(colony)
    obs = Observable(name='length', raw='Long axis (L)')
    colony.decompose(seed=357)  # so that decomposition is always identical
    myplot = SamplePlot([colony, ], parser=parser)
    myplot.make_plot(obs)
    myplot.save(user_bname='sample', add_obs=False, extension='.png')

if __name__ == '__main__':
    location = '~/Boulot/Temp/seggerdir'  # this is a supersegger root folder
    exp = Experiment(path=location, filetype='supersegger')
    cell = iterate(exp)
#
#    exp._count_items()

    plot_samples(exp)