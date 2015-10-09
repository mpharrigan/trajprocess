import os

from nose import with_setup
import mdtraj

from .mock1 import generate_project, cleanup, generate_bw
from .mock2 import mock_project as mock2, cleanup as cleanup2
from trajprocess.files import Project, record
from trajprocess import process

import numpy as np
import json

import subprocess


@with_setup(generate_project, cleanup)
def test_cnv():
    project = Project("p1234", "data/PROJ1234", 'xa4')
    raw_infos = list(project.get_infos())
    nfo_infos = list(map(record(project.nfo), raw_infos))
    cat_infos = list(map(record(project.cat), nfo_infos))
    cnv_infos = list(map(record(project.cnv), cat_infos))

    assert len(cnv_infos) > 0

    for info in cnv_infos:
        assert os.path.exists(info['cnv']['xtc_out'])
        assert os.path.exists(info['cnv']['nc_out'])

        with mdtraj.open(info['cnv']['xtc_out']) as tfile:
            with mdtraj.open(info['cnv']['nc_out']) as dcdfile:
                xyz, time, step, box = tfile.read()
                print("Shape", xyz.shape)
                assert xyz.shape == (21, 22, 3), xyz.shape

                xyz_nc, time, lengths, angles = dcdfile.read()
                np.testing.assert_array_equal(xyz * 10, xyz_nc)


@with_setup(generate_bw, cleanup)
def test_cnv_bw():
    project = Project("v1", "data/v1", 'bw')
    raw_infos = list(project.get_infos())
    nfo_infos = list(map(record(project.nfo), raw_infos))
    cat_infos = list(map(record(project.cat), nfo_infos))
    cnv_infos = list(map(record(project.cnv), cat_infos))

    assert len(cnv_infos) > 0

    for info in cnv_infos:
        assert os.path.exists(info['cnv']['xtc_out'])
        assert "cnv.xtc" in info['cnv']['xtc_out'], info['cnv']['xtc_out']

        with mdtraj.open(info['cnv']['xtc_out']) as tfile:
            xyz, time, step, box = tfile.read()
            print("Shape", xyz.shape)
            assert xyz.shape == (11, 22, 3), xyz.shape


@with_setup(mock2, cleanup2)
def test_cnv_nc():
    os.remove("processed/p9761/24/7/cnv.nc")
    with open("processed/p9761/24/7/info.json") as f:
        info = json.load(f)
    process.cnv_to_nc(info, chunk=2)

    top = 'tops-p9712/{top[struct]}.prmtop'.format(**info)
    trj1 = mdtraj.load(info['cnv']['xtc_out'], top=top)
    trj2 = mdtraj.load(info['cnv']['nc_out'], top=top)

    np.testing.assert_array_almost_equal(trj1.xyz, trj2.xyz)
    np.testing.assert_array_almost_equal(trj1.unitcell_vectors,
                                         trj2.unitcell_vectors)
    np.testing.assert_array_almost_equal(trj1.time, trj2.time)


@with_setup(mock2, cleanup2)
def test_nc_cpptraj():
    os.remove("processed/p9761/24/7/cnv.nc")
    with open("processed/p9761/24/7/info.json") as f:
        info = json.load(f)
    process.cnv_to_nc(info, chunk=2)

    top = 'tops-p9712/{top[struct]}.prmtop'.format(**info)
    out = "{workdir}/test.nc".format(**info['path'])
    subprocess.check_call([
        'cpptraj',
        '-p', top,
        '-y', info['cnv']['nc_out'],
        '-x', out
    ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    trj1 = mdtraj.load(out, top=top)
    trj2 = mdtraj.load(info['cnv']['nc_out'], top=top)

    np.testing.assert_array_almost_equal(trj1.xyz, trj2.xyz)
    np.testing.assert_array_almost_equal(trj1.unitcell_vectors,
                                         trj2.unitcell_vectors)
    #np.testing.assert_array_almost_equal(trj1.time, trj2.time)
