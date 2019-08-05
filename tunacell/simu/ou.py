#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Classes and functions to simulate an Ornstein-Uhlenbeck process for the instantaneous growth rate


See also
--------
:mod:`tunacell.simu.base` for base classes to parameterize numerical simulations

"""
from __future__ import print_function
from builtins import input  # future package

import os
import shutil
import logging
import datetime
import uuid

from tqdm import tqdm

import numpy as np

from tunacell.base.experiment import Experiment
from tunacell.base.container import Container
from tunacell.base.colony import Colony

from tunacell.io.metadata import Metadata

from tunacell.simu.base import Ecoli, SimuParams, DivisionParams, SampleInitialSize


def run_ou_simulation(
    simuParams, divParams, bsParams, ouParams, where, label, force=False
):
    """API function to run and export OU numerical simulation

    Parameters
    ----------
    simuParams : :class:`tunacell.simu.base.SimuParams` instance
        sets general simulation parameters
    divParams : :class:`tunacell.simu.base.DivisionParams` instance
        sets cell division process
    bsParams : :class:`tunacell.simu.base.SampleInitialSize` instance
        sets how initial cell size is sampled in simulation
    ouParams : :class:`OUParams`
        sets Ornstein-Uhlenbeck parameters for simulation of growth rate
    where : str
        parent directory where the experiment will be exported
    label : str
        name of the experiment to be exported
    force : bool (default False)
        whether to overide previous same-label experiment folder
    """
    exp = OUSimulation(
        label=label,
        simuParams=simuParams,
        divisionParams=divParams,
        ouParams=ouParams,
        birthsizeParams=bsParams,
    )

    # check that experiment has been saved before
    ans = "go"
    current_name = label
    exp_path = os.path.join(where, label)
    if force:
        if os.path.exists(exp_path):
            shutil.rmtree(exp_path)
    while os.path.exists(exp_path) and ans != "a":
        print("Experiment {} already exists.".format(current_name))
        ans = input("Override [o], Change experiment name [c], Abort [a]: ")
        if ans == "c":
            new_name = current_name
            while new_name == current_name:
                new_name = input("NEW NAME: ")
            exp.label = new_name
            exp_path = os.path.join(where, new_name)
        # when overriding, need to erase everything first
        elif ans == "o":
            shutil.rmtree(exp_path)

    # export except if process has been aborted
    if ans != "a":
        exp.raw_text_export(path=where)


class OUSimulation(Experiment):
    """Equivalent of Experiment class for simulation.

    An instance of OUSimulation contains all information to simulate an
    Ornstain-Uhlenbeck process in dividing cells.

    Parameters
    ----------
    simuParams : SimuParams instance
        provides simulation parameters (number of containers, number of
        colonies per container, timing parameters)
    divisionParams : DivisionParams instance
        sets the division timing process
    ouParams : OUParams instance
        set the parameters for the Ornstein Uhlenbeck process
    """

    def __init__(
        self,
        label=None,
        simuParams=None,
        divisionParams=None,
        ouParams=None,
        birthsizeParams=None,
        filter_set=None,
    ):
        today = datetime.datetime.today()
        self.date = today
        if simuParams is None:
            simuParams = SimuParams()
            print("Using default SimuParams:\n{}".format(simuParams))
        if divisionParams is None:
            divisionParams = DivisionParams()
            print("Using default DivisionParams:\n{}".format(divisionParams))
        if ouParams is None:
            ouParams = OUParams()
            print("Using default OUParams:\n{}".format(ouParams))
        if birthsizeParams is None:
            birthsizeParams = SampleInitialSize()
        if label is None:
            self._label = "simu_{}".format(today.strftime("%Y-%m-%d_%H-%M-%S"))
        else:
            self._label = label
        self.abspath = "{}".format(hex(id(self)))  # Experiment compatibility
        self.filetype = "simu"
        self.datatype = [
            ("time", "f8"),
            ("ou", "f8"),
            ("ou_int", "f8"),
            ("exp_ou_int", "f8"),
            ("cellID", "u2"),
            ("parentID", "u2"),
        ]

        self.containers = []  # there are no file
        self.simuParams = simuParams
        self.divisionParams = divisionParams
        self.ouParams = ouParams
        self.birthsizeParams = birthsizeParams

        self._set_metadata()
        # set filterset
        self.fset = filter_set
        return

    def _set_metadata(self):

        dic = {
            "experiment": self.label,  # mandatory
            "period": self.simuParams.period,  # mandatory
            "date": self.date.strftime("%Y-%m-%d"),  # optional
        }

        # optional, but convenient: paramaters as sub-dicts
        dic["simu_params"] = {k: v for k, v in self.simuParams.content}
        dic["division_params"] = {k: v for k, v in self.divisionParams.content}
        dic["ornstein_uhlenbeck_params"] = {k: v for k, v in self.ouParams.content}
        dic["birth_size_params"] = {k: v for k, v in self.birthsizeParams.content}

        #        dics = {key: {self.label: value} for key, value in content}
        #        self.metadata = pd.DataFrame(dics)

        self.metadata = Metadata([dic])

        return

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, lab):
        if not isinstance(lab, str):
            lab = repr(lab)
        self._label = lab
        # reset metadata in case label is changed
        self._set_metadata()
        return

    def info(self):
        msg = (
            "Ornstein Uhlenbeck simulation on dividing cells" + "\n"
            "-----------------------------------------------" + "\n"
            "{}".format(self.simuParams) + "\n"
            "{}".format(self.divisionParams) + "\n"
            "{}".format(self.ouParams)
        )
        return msg

    def iter_containers(
        self,
        size=None,
        read=True,
        build=True,  # useless options
        prefilt=None,  # only used for compatibility
        extend_observables=True,  # idem
        report_NaNs=True,  # idem
        shuffle=False,
    ):  # idem
        if size is None:
            size = self.simuParams.nbr_container
        # determine number of digits
        ndigits = len("{:d}".format(size))
        # random seed
        np.random.seed(self.simuParams.seed)
        # start loop over containers
        for index in tqdm(range(size), desc="simulation"):
            label = "container_{nbr:0{digits:d}d}".format(digits=ndigits, nbr=index + 1)
            yield OUContainer(
                self,
                label=label,
                simuParams=self.simuParams,
                divisionParams=self.divisionParams,
                ouParams=self.ouParams,
                birthsizeParams=self.birthsizeParams,
            )
        return


class OUContainer(Container):
    """Ornstein-Uhlenbeck simulation container.

    This is subclassed from :class:`tunacell.base.container.Container`,
    only ``__init__`` changes.


    Parameters
    ----------
    simu : OUSimulation instance
    label : str
        label for the container
    simuParams : SimuParams instance
    divisionParams : DivisionParams instance
    ouParams : OUParams instance
    birthsizeParams : SampleInitialSize instance
    """

    def __init__(
        self,
        simu,
        label=None,
        simuParams=None,
        divisionParams=None,
        ouParams=None,
        birthsizeParams=None,
    ):
        """Runs the simulation upon call.

        """
        self.exp = simu  # Container compatibility
        self.abspath = "{}".format(hex(id(self)))  # Container compatibility
        self.filetype = "simulations"
        self.datatype = simu.datatype
        self.metadata = simu.metadata

        if label is not None:
            self.label = label
        else:
            self.label = str(uuid.uuid1())[:8]
        self.cells = []
        self.trees = []
        nodes = []
        trees = []
        count = 0
        for samp in range(simuParams.nbr_colony_per_container):
            colony, count = ou_tree(
                ouParams,
                divisionParams,
                birthsizeParams,
                count=count + 1,
                tstart=simuParams.start,
                tstop=simuParams.stop,
                dt=simuParams.period,
            )
            colony.container = self
            trees.append(colony)
            for node in colony.all_nodes():
                nodes.append(node)
        self.trees = trees
        self.cells = nodes
        return


class OUParams(object):
    """Class to store Ornstein-Uhlenbeck parameters.

    Parameters
    ----------
    target : float (default np.log(2)/60)
        target value for growth rate, default value is one doubling time per hour in units of min^-1
    spring : float (default 1/30)
        spring constant (inverse of autocorrelation time), default value is 1/30 per minute
    noise : float (default 2.0/30.0 * (np.log(2.) / 600.)**2)
        sets the noise intensity, default value is chosen to get a 10% standard deviation at steady state
        with the two other default parameters

    See also
    --------
    Gillespie, D.T., Phys Rev E, vol 54, pp 2084-2091 (1996)
    """

    def __init__(
        self,
        target=np.log(2.0) / 60,
        spring=1.0 / 30.0,
        noise=2.0 / 30.0 * (np.log(2.0) / 600.0) ** 2,
    ):
        self.target = target
        self.spring = spring
        self.noise = noise
        # metadata helper
        self.content = [
            ("target", float(target)),
            ("spring", float(spring)),
            ("noise", float(noise)),
        ]
        return

    def __repr__(self):
        lab = ""
        lab += "Ornstein-Uhlenbeck parameters:\n"
        lab += "Equation: dx/dt = -k * (x - x_0) + xi(t)\n"
        lab += "with <xi(t)xi(s)> = eta^2 delta(t-s)\n"
        lab += "* x_0, target value: {}\n".format(self.target)
        lab += "* k, spring constant: {}\n".format(self.spring)
        lab += "* eta^2, noise intensity: {}\n".format(self.noise)
        return lab


class OUsteps(object):
    """Class that computes Gillespie-like prefactors for OU updates.

    Parameters
    ----------
    params : OUParams instance
        set of parameters for OU process
    dt : float
        time interval for update.

    Attributes
    ----------
    mu : float
        target value
    sigma : float
    sigma_y : float
    kappa : float

    See also
    --------
    Gillespie, D.T., Phys Rev E, vol 54, pp 2084-2091 (1996)
    """

    def __init__(self, params, dt=1.0):
        # Gillespie-like parameters
        k = params.spring
        c = params.noise
        mu = np.exp(-k * dt)
        self.mu = mu
        sx2 = c / (2.0 * k) * (1.0 - np.exp(-2.0 * k * dt))
        sy2 = c / (k ** 3) * (k * dt - 2.0 * (1.0 - mu) + 0.5 * (1.0 - mu ** 2))
        kappa = c * (1.0 - mu) ** 2 / (2.0 * k ** 2)
        self.sigma = np.sqrt(sx2)
        self.sigma_y = np.sqrt(sy2)
        self.kappa = kappa
        return


def ou_track(params, dt=1.0, steps=10, start=0.0, y_start=1.0):
    """Sample OU process X(t) and its integral Y(t).

    Beware that X(t) is the zero mean Ornstein-Uhlenbeck process.

    Parameters
    ----------
    params : OUParams instance
    dt : float
        time interval between two steps
    steps : int
        number of steps in track
    start : float
        value of OU process at initial step (step 0)
    y_start : float
        value of integrated OU process at initial step (step 0)

    Returns
    -------
    (OU X(t), OU Y(t))
    OU X(t) : ndarray
        Ornstein-Uhlenbeck process values sampled at time interval dt
    OU Y(t) : ndarray
        Ornstein-Uhlenbeck integrated process values

    Notes
    -----
    X(t) is defined by:

    ..math::
        \frac{d X(t)}{dt} = - k * X(t) + c^{1/2} \Gamma(t)

    where :math:`Gamma(t)` is a Gaussian white noise with zero mean and
    covariance
    ..math::
        \langle \Gamma(t) \Gamma(s) \rangle = \delta(t-s)

    where :math:`\delta(t-s)` is the Dirac delta function.

    The integrated process Y(t) follows

    ..math::
        \frac{d Y(t)}{dt} = X(t)

    See also
    --------
    Gillespie, D.T., Phys Rev E, vol 54, pp 2084-2091 (1996)
    """

    xs = np.zeros(steps + 1, dtype="f8")
    ys = np.zeros(steps + 1, dtype="f8")

    ns = np.random.normal(0.0, 1.0, size=(steps, 2))
    xs[0] = start
    ys[0] = y_start

    k = params.spring  # spring constant

    ou = OUsteps(params, dt=dt)  # calling to get Gillespie-like step params
    mu = ou.mu

    sigma_x = ou.sigma
    sigma_y = ou.sigma_y
    kappa = ou.kappa

    for i, (n1, n2) in enumerate(ns, start=1):
        # Gillespie update
        xs[i] = xs[i - 1] * mu + sigma_x * n1
        ys[i] = ys[i - 1] + (
            xs[i - 1] * (1.0 - mu) / k
            + np.sqrt((sigma_y) ** 2 - (kappa ** 2 / sigma_x ** 2)) * n2
            + kappa / sigma_x * n1
        )
    return xs, ys


def ou_tree(
    ouparams, divparams, birthsizeparams, count=None, tstart=0.0, tstop=300.0, dt=5.0
):
    """Generates recursively OU process on dividing cells.

    Parameters
    ----------
    ouparams : OUParams instance
    divparams : DivisionParams instance
    birthsizeparams : SampleInitialSize instance
    count : int
        label for root cell
    tstart : float
        time at which simulation starts
    tstop : float
        time at which simulation stops
    dt : float
        time interval at which value of OU process are recorded. Such value is
        used to reject smaller interdivision times.

    Returns
    -------
    tree : Colony instance
        in which process is stored in each node Ecoli.data
    count : int
        last cell label used
    """
    if count is not None:
        rootid = str(count)  # labeling by integers (exported as strings)
    else:
        rootid = None  # automatic labeling
    root = root_cell(
        ouparams, divparams, birthsizeparams, identifier=rootid, tstart=tstart, dt=dt
    )
    tree = Colony()
    tree.add_node(root)
    count = add_recursive_branch(
        root,
        tree,
        count=count,
        tstart=tstart,
        tstop=tstop,
        dt=dt,
        ouparams=ouparams,
        divparams=divparams,
    )

    return tree, count


def root_cell(
    ouparams, divparams, birthsizeparams, identifier=None, tstart=0.0, dt=5.0
):
    """Set state for root cell, that initialize a sample.

    Parameters
    ----------
    ouparams : OUParams instance
        store OU parameters: target, spring, noise
    divparams : DivisionParams instance
        store information to generate random cell cycle duration
    birthsizeParams : SampleInitialSize instance
        sets sampling of root cell birth size
    identifier : str or int (default None)
        identifier to give to root cell. Give an integer to increase +1.
    tstart : float (default 0.)
        start of recording simulations
    dt : float (default 5.)
        time between two successive acquisitions; smaller interdivision times
        are rejected.
    birth_size : float (default 1.)
        size at birth for the root cell
    """
    # start with OU equilibrium sample for alpha
    equilibrium_mean = ouparams.target
    equilibrium_std = np.sqrt(ouparams.noise / (2.0 * ouparams.spring))
    birth_growth_rate = -1
    while birth_growth_rate < 0:  # prevent surprises...
        birth_growth_rate = np.random.normal(
            loc=equilibrium_mean, scale=equilibrium_std
        )
    birth_size = birthsizeparams.rv()

    if divparams.use_growth_rate == "parameter":
        growth_rate = ouparams.target
    elif divparams.use_growth_rate == "birth":
        growth_rate = birth_growth_rate
    # reject lifetime smaller than dt
    counter = 0
    lifetime_root = -1
    while lifetime_root < dt and counter < 1000:
        lifetime_root = divparams.rv(birth_size, growth_rate)
        if lifetime_root < dt:
            logging.info(
                "Rejecting interdivision time {} < period {}".format(lifetime_root, dt)
            )
        counter += 1
    if counter == 1000:
        raise ValueError("Interdivision time always larger than acquisition period.")

    age_root = np.random.uniform()
    birth_root = tstart - age_root * lifetime_root

    root = Ecoli(
        identifier=identifier,
        parent=None,
        birth_time=birth_root,
        lifetime=lifetime_root,
    )

    root.birth_value = (birth_growth_rate, np.log(birth_size))

    return root


def add_recursive_branch(
    ecoli,
    tree,
    count=None,
    tstart=0.0,
    tstop=300.0,
    dt=5.0,
    ouparams=OUParams(),
    divparams=DivisionParams(),
):
    """Main function for generating the tree and simulated process.
    """
    x_start, y_start = ecoli.birth_value
    t_birth = ecoli.birth_time
    t_div = ecoli.division_time
    if count is None:  # in this case, identifier is a 36 char string from uuid
        ecoli.tag = ecoli.identifier[:8]
        if isinstance(ecoli.identifier, str):
            idtype = "U36"  # unicode string
        else:
            idtype = "S36"  # byte string

        def store_id(identifier):
            return identifier

    else:
        ecoli.tag = ecoli.identifier
        idtype = "u2"

        # in this case, identifier is a string made from integer
        def store_id(identifier):
            return int(identifier)

    #    print 'add {}, birth {:.1f}, div {:.1f}'.format(ecoli.tag,
    #                                                    ecoli.birth_time,
    #                                                    ecoli.division_time)
    # print(idtype)
    first_rec_time = tstart + (np.floor((t_birth - tstart) / dt) + 1.0) * dt
    first_rec_time = max(tstart, first_rec_time)  # selection for root cell

    rec_times = np.arange(first_rec_time, min(t_div, tstop), dt)
    # (this was a selection for leaves: the ones that cross tstop)

    cid = store_id(ecoli.identifier)
    if ecoli.bpointer is None:
        pid = store_id("0")
    else:
        pid = store_id(ecoli.bpointer)

    ecoli.rec_times = rec_times
    # Check that there is at least one recording time
    if len(rec_times):
        # first step from birth to first recording time
        (xs, ys) = ou_track(
            ouparams,
            dt=rec_times[0] - t_birth,
            steps=1,
            start=x_start - ouparams.target,
            y_start=y_start,
        )
        first_x = xs[-1]
        first_y = ys[-1]
        # record values
        (xs, ys) = ou_track(
            ouparams, dt=dt, steps=len(rec_times) - 1, start=first_x, y_start=first_y
        )
        x_rec_values = xs + ouparams.target
        y_rec_values = ys + ouparams.target * (rec_times - t_birth)
        length_like_values = np.exp(y_rec_values)
        last_x = xs[-1]
        last_y = ys[-1]
        rec_ids = len(rec_times) * [cid]
        rec_pids = len(rec_times) * [pid]
        ecoli.data = np.zeros(
            len(rec_times),
            dtype=[
                ("time", "f8"),
                ("ou", "f8"),
                ("ou_int", "f8"),
                ("exp_ou_int", "f8"),
                ("cellID", idtype),
                ("parentID", idtype),
            ],
        )
        ecoli.data["time"] = rec_times
        ecoli.data["ou"] = x_rec_values
        ecoli.data["ou_int"] = y_rec_values
        ecoli.data["exp_ou_int"] = length_like_values
        ecoli.data["cellID"] = rec_ids
        ecoli.data["parentID"] = rec_pids

        last_dt = t_div - rec_times[-1]
    # otherwise, data is set to empty array, but we update cycle bounds
    else:
        last_x = x_start - ouparams.target
        last_y = y_start
        last_dt = t_div - t_birth
        # empty array since no recording time between birth and division
        ecoli.data = np.array(
            [],
            dtype=[
                ("time", "f8"),
                ("ou", "f8"),
                ("ou_int", "f8"),
                ("exp_ou_int", "f8"),
                ("cellID", idtype),
                ("parentID", idtype),
            ],
        )

    # last step
    (xs, ys) = ou_track(ouparams, dt=last_dt, steps=1, start=last_x, y_start=last_y)
    ecoli.division_value = (
        xs[-1] + ouparams.target,
        ys[-1] + ouparams.target * (t_div - t_birth),
    )

    # create two daughter cells if time has not reached tmax
    if t_div < tstop:
        new_birth_size = np.exp(ecoli.division_value[1]) / 2.0  # symmetric division
        new_alpha = ecoli.division_value[0]
        for i in range(2):
            if divparams.use_growth_rate == "parameter":
                alpha = ouparams.target
            elif divparams.use_growth_rate == "birth":
                alpha = new_alpha
            # compute cell interdivision time and accept if larger than period
            counter = 0
            lt = -1
            while lt < dt and counter < 1000:
                lt = divparams.rv(new_birth_size, alpha)
                if lt < dt:
                    logging.info(
                        "Rejecting interdivision time {} < period {}".format(lt, dt)
                    )
                counter += 1
            if counter == 1000:  # stupid control
                raise ValueError(
                    "Interdivision time always smaller than acquisition period"
                )
            # moving on
            if count is not None:
                count += 1
                newid = str(count)
            else:
                newid = None  # will create identifier as uuid item
            necoli = Ecoli(identifier=newid, parent=ecoli, lifetime=lt)
            # update necoli.birth_value for y process
            necoli.birth_value = (new_alpha, np.log(new_birth_size))
            tree.add_node(necoli, parent=ecoli.identifier)
            count = add_recursive_branch(
                necoli,
                tree,
                count=count,
                tstart=tstart,
                tstop=tstop,
                dt=dt,
                ouparams=ouparams,
                divparams=divparams,
            )

    return count  # this is the last taken integer + 1
