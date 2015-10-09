import os
import subprocess
import tarfile
import json

from nose import with_setup
import mdtraj

from trajprocess.postprocess import _norm_cpptraj
import trajprocess.postprocess
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

    assert os.path.exists(info['cnv']['xtc_out'])


@with_setup(mock_project, cleanup)
def test_stp_trek():
    with open("processed/p9761/24/7/info.json") as f:
        info = json.load(f)
    info['cnv']['nc_out'] = "{workdir}/cnv.nc".format(**info['path'])

    # Make sure we can overwrite this
    with open("{workdir}/cpptraj.tar.gz".format(**info['path']), 'wb') as f:
        f.write(b'WHATUP')

    info = trajprocess.postprocess.stp_trek(info)

    assert tarfile.is_tarfile('{workdir}/cpptraj.tar.gz'.format(**info['path']))
    assert not os.path.exists('{workdir}/cpptraj'.format(**info['path']))

    traj = mdtraj.load(info['stp']['nc_out'], top=info['stp']['prmtop'])
    assert traj.n_atoms == 30962
    assert len(traj) == 7
