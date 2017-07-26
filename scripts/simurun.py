#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Script to run numerical simulations.
"""
from __future__ import print_function

import argparse
import os
import shutil
import numpy as np
from tuna.simu.main import SimuParams, DivisionParams
from tuna.simu.ou import OUParams, OUSimulation

# Arguments
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--path', type=str,
                    help='Parent directory in which simulation is stored',
                    default='~/tmptuna')
parser.add_argument('-l', '--label', type=str,
                    help='Label of the experiment/simulation',
                    default='simutest')
parser.add_argument('-s', '--samples', type=int,
                    help='Number of simulated container samples',
                    default=100)
parser.add_argument('--start', type=float,
                    help='Time at which simulation starts',
                    default=0.)
parser.add_argument('--stop', type=float,
                    help='Time at which simulation stops',
                    default=180.)
parser.add_argument('--period', type=float,
                    help=('Time period between two consecutive time-lapse'
                          ' acquisitions'),
                    default=5.)
args = parser.parse_args()


path = os.path.abspath(os.path.expanduser(args.path))
if not os.path.exists(path):
    os.makedirs(path)
print('Path: {}'.format(path))
print('Label: {}'.format(args.label))
# %% SET PARAMETERS FOR NUMERICAL SIMULATION

# general simulation parameters
tstart = args.start  # initial time of simulation
tstop = args.stop  # final time of simulation (data is not recorded beyond this)
dt = args.period  # time interval for data recording (sampling of process)
nsamples = args.samples  # number of containers that will be created by simulation

tree_per_cont = 2  # number of trees per container

# OU parameters
target_value = np.log(2.)/60.  # average value of the random process
spring = 1./30.  # spring constant (brings back random walker to target value)
noise_intensity = 2. * spring * (target_value/10.)**2   # noise intensity for random walk

# DIVISION TIMING
interdivision_mean = 60.  # average of interdivision time
interdivision_std = 6.  # standard deviation for interdivision time
minimum = dt  # minimum value of interdivision time should be >= dt

# %% Instantiating parameters objects

simuParams = SimuParams(nbr_container=nsamples,
                        nbr_colony_per_container=tree_per_cont,
                        start=tstart, stop=tstop, interval=dt)

ouParams = OUParams(target=target_value, spring=spring,
                    noise=noise_intensity)

divParams = DivisionParams(mean=interdivision_mean, std=interdivision_std)

# %% SIMULATION: initiate, define, run/write
np.random.seed()  # seed for the random number generator

exp = OUSimulation(label=args.label,
                   simuParams=simuParams, divisionParams=divParams,
                   ouParams=ouParams)

# check that experiment has been saved before
ans = 'go'
current_name = args.label
exp_path = os.path.join(path, args.label)
while os.path.exists(exp_path) and ans != 'a':
    print('Experiment {} already exists.'.format(current_name))
    ans = raw_input('Override [o], Change experiment name [c], Abort [a]: ')
    if ans == 'c':
        new_name = current_name
        while new_name == current_name:
            new_name = raw_input('NEW NAME: ')
        exp.label = new_name
        exp_path = os.path.join(path, new_name)
    # when overriding, need to erase everything first
    elif ans == 'o':
        shutil.rmtree(exp_path)

# export except if process has been aborted
if ans != 'a':
    exp.raw_text_export(path=path)
