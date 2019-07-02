#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
This module indirectly defines the command-line: tunasimu (see entry point in setup.py)

It executes numerical simulations from command-line. Usage:

    tunasimu [[--dir <path>]
               [-l|--label <exp-name>]
               [-s|--seed <seed>]
               [-f|--force]
               [-s|--samples <number-of-samples>]
               [--alpha <growth-rate-target>]
               [--alpha_sd <standard-deviation-for-alpha]
               [--alpha_autocorr <autocorrelation time for alpha>]
               [--div_lambda <lambda division parameter>]
               [--div_size_cutoff <target size for division>]
               [--div_use_alpha {'parameter', 'birth'}]
               [--div_fluctuation_mode {'gamma', 'none'}]
               [--div_relative_fluctuations <ratio sd/mean>]
               [--initial_birth_size_mode {'fixed', 'lognormal'}]
               [--initial_birth_size_sigma <sigma for initial birth size>]
               [--start <starting-time-value>]
               [--stop <stoping-time-value>]
               [--period <time-interval-between-acquisitions>]
              ]

All parameters have default values (see code).

(the execution can be stated as:

    python -m tunacell.simu.run [options...]
"""
from __future__ import print_function

import argparse
import os
import yaml

import numpy as np

from tunacell.simu.ou import (OUParams,
                              SimuParams, DivisionParams, SampleInitialSize,
                              run_ou_simulation)



def setup_args():
    """Configure argparse parser to get all parameters for numerical simulation
    """
    # Arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # files
    parser.add_argument('--dir', type=str,
                        help='Parent directory in which simulation is stored',
                        default='~/tmptunacell')
    parser.add_argument('-l', '--label', type=str,
                        help='Label of the experiment/simulation',
                        default='simutest')
    parser.add_argument('-c', '--conf', type=str,
                        help='Configuration file (e.g. "metadata.yml")',
                        default=None)
    parser.add_argument('-s', '--seed', type=int,
                        help='Random generator seed (use a specific number for reproducibility)',
                        default=None)
    parser.add_argument('-f', '--force',
                        help='Force writting new data (erase old simulations)',
                        action='store_true')
    parser.add_argument('-N', '--N_containers', type=int,
                        help='Number of simulated container samples',
                        default=100)
    parser.add_argument('-t', '--tree_per_container', type=int,
                        help='Number of trees per container',
                        default=2)
    # OU parameters
    parser.add_argument('--alpha', type=float,
                        help='Growth rate target (in [/min])',
                        default=np.log(2.)/60.)  # corresponding to size doubling time = 1 hour, in /min
    parser.add_argument('--alpha_sd', type=float,
                        help='Equilibrium standard deviation for alpha (in [/min])',
                        default=np.log(2.)/(10.*60.))  #1/10 of default target, in /min
    parser.add_argument('--alpha_autocorr', type=float,
                        help='Autocorrelation time for alpha (in [min]',
                        default=30.)  # in mins
    # division parameters
    parser.add_argument('--div_lambda', type=float,
                        help='Lambda parameter for interdivision timing (between 0 and 1)',
                        default=0)  # timer
    parser.add_argument('--div_size_cutoff', type=float,
                        help='Division size cutoff (in arbitrary units)',
                        default=2.)
    parser.add_argument('--div_use_alpha', type=str,
                        help='Use parameter, or birth value growth rate',
                        default='parameter')
    parser.add_argument('--div_fluctuation_mode', type=str,
                        help='[gamma, none] whether to add or not fluctuations to interdivision time',
                        default='gamma')
    parser.add_argument('--div_relative_fluctuations', type=float,
                        help='Ratio sd/predicted linear interdivision time',
                        default=.1)
    
    # initial birth size sampling
    parser.add_argument('--birth_size_mean', type=float,
                        help='Initial mean birth size', default=1.)
    parser.add_argument('--birth_size_mode', type=str,
                        help='[fixed, lognormal]', default='fixed')
    parser.add_argument('--birth_size_relative_fluctuations', type=float,
                        help='Relative fluctuations in the lognormal mode',
                        default=0.1)
    
    # simulation parameters
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
#    print(vars(args).keys())
    return args


def main():
    """Setup numerical simulation parameters from command line arguments
    
    Parameters
    ----------
    args : :class:`argparse.Namespace` instance
        any Python object is valid as it contains following attributes:
        ['dir', 'label', 'random_seed', 'force', 'samples',
        'alpha', 'alpha_sd', 'alpha_autocorr', 'div_lambda', 'div_size_cutoff',
        'div_use_alpha', 'div_fluctuation_mode', 'div_relative_fluctuations',
        'initial_birth_size_mode', 'initial_birth_size_sigma', 'start', 'stop',
        'period']
    """
    args = setup_args()

    path = os.path.abspath(os.path.expanduser(args.dir))
    if not os.path.exists(path):
        os.makedirs(path)
    print('Path: {}'.format(path))
    print('Label: {}'.format(args.label))
    # SET PARAMETERS FOR NUMERICAL SIMULATION
    
    # check that a configuration file is provided
    if args.conf is not None:
        with open(args.conf, 'r') as f:
            y = yaml.load(f)
        simuParams = SimuParams(**y['simu_params'])
        ouParams = OUParams(**y['ornstein_uhlenbeck_params'])
        divParams = DivisionParams(**y['division_params'])
        birthsizeParams = SampleInitialSize(**y['birth_size_params'])
    # otherwise read default
    else:
        # general simulation parameters
        tstart = args.start  # initial time of simulation
        tstop = args.stop  # final time of simulation (data is not recorded beyond this)
        dt = args.period  # time interval for data recording (sampling of process)
        ncontainers = args.N_containers  # number of containers that will be created by simulation
        ntrees = args.tree_per_container  # number of trees per container
        seed = args.seed
        
        simuParams = SimuParams(nbr_container=ncontainers,
                                nbr_colony_per_container=ntrees,
                                start=tstart, stop=tstop, period=dt,
                                seed=seed)
        
        # OU parameters
        target_value = args.alpha  # average value of the random process
        spring = 1./args.alpha_autocorr  # spring constant (brings back random walker to target value)
        noise_intensity = 2. * spring * (args.alpha_sd)**2   # noise intensity for random walk
        
        ouParams = OUParams(target=target_value, spring=spring, noise=noise_intensity)
        
        # DIVISION TIMING
        divParams = DivisionParams(div_lambda=args.div_lambda,
                                   div_size_target=args.div_size_cutoff,
                                   div_mode=args.div_fluctuation_mode,
                                   div_sd_to_mean=args.div_relative_fluctuations,
                                   use_growth_rate=args.div_use_alpha)
        
        # INITIAL BIRTH SIZE
        birthsizeParams = SampleInitialSize(birth_size_mean=args.birth_size_mean,
                                            birth_size_mode=args.birth_size_mode,
                                            birth_size_sd_to_mean=args.birth_size_relative_fluctuations)
    
    run_ou_simulation(simuParams, divParams, birthsizeParams, ouParams,
                      path, args.label, force=args.force)
    
    
if __name__ == '__main__':
    main()
