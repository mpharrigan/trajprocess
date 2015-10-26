import os
import subprocess
import tarfile
import json
import numpy as np

from nose import with_setup
import mdtraj

from trajprocess.postprocess import _norm_cpptraj
import trajprocess.postprocess
from trajprocess.files import Postprocess, process_post
from .mock2 import mock_project, cleanup


def test_norm_cpptraj():
    assert _norm_cpptraj(":WAT") == 'mol-wat'
    assert _norm_cpptraj(":MY") == 'mol-my'
    assert _norm_cpptraj("@Na+") == "atm-na-pl"
    assert _norm_cpptraj("@Cl-") == "atm-cl-min"


def test_cpptraj_exists():
    subprocess.check_call(['cpptraj', '--version'], stdout=subprocess.DEVNULL)


@with_setup(mock_project, cleanup)
def test_mock2():
    assert os.path.exists("processed")
    with open("processed/p9761/24/7/info.json") as f:
        info = json.load(f)

    assert os.path.exists(info['cnv2']['gens'][0])


@with_setup(mock_project, cleanup)
def test_trek():
    # setup
    with open("processed/p9761/24/7/info.json") as f:
        info = json.load(f)

    info = trajprocess.postprocess.stp(info, 'trek')

    # check stp cleanup
    assert not os.path.exists('{workdir}/stp/0/'.format(**info['path']))

    # check stp results
    traj = mdtraj.load(info['stp']['gens'][0], top=info['stp']['outtop'])
    assert traj.n_atoms == 30962
    assert len(traj) == 7

    # do ctr
    info = trajprocess.postprocess.ctr(info, "trek")

    # check ctr info
    assert not os.path.exists("{workdir}/cpptraj.tmp".format(**info['path']))
    assert not os.path.exists(
        "{workdir}/ctr/cpptraj.tmp".format(**info['path']))
    traj2 = mdtraj.load(info['ctr']['gens'][0], top=info['stp']['outtop'])

    # check ctr results
    # Trek has 518 protein residues
    pairs = np.random.randint(0, 518, (20, 2))
    cont1, _ = mdtraj.compute_contacts(traj, pairs)
    cont2, _ = mdtraj.compute_contacts(traj2, pairs)

    np.testing.assert_array_almost_equal(cont1, cont2, decimal=4)


@with_setup(mock_project, cleanup)
def test_process_post():
    # setup
    with open("processed/p9761/24/7/info.json") as f:
        info = json.load(f)
    infos = [info]

    infos = process_post(Postprocess('trek'), infos)
    info, *_ = infos
    traj = mdtraj.load(info['ctr']['gens'][0], top=info['stp']['outtop'])
    assert traj.n_atoms == 30962

    with open("processed/p9761/24/7/info.json") as f:
        disk_info = json.load(f)

    assert disk_info == info
