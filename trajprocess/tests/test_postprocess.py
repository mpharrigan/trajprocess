from trajprocess.postprocess import _norm_cpptraj
import trajprocess.postprocess
import os

from nose import with_setup
import mdtraj
from trajprocess.files import Project, record
import subprocess

from .mock2 import mock_project, cleanup

import json


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
def test_stp():
    with open("processed/p9761/24/7/info.json") as f:
        info = json.load(f)
    info['cnv']['nc_out'] = "{workdir}/cnv.nc".format(**info['path'])

    trajprocess.postprocess.stp_trek(info)
