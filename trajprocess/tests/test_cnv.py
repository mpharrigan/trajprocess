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
    cnv1_infos = list(map(record(project.cnv1), nfo_infos))
    cnv2_infos = list(map(record(project.cnv2), cnv1_infos))

    assert len(cnv2_infos) > 0

    for info in cnv2_infos:
        for g in range(len(info['cnv2']['gens'])):
            assert os.path.exists(info['cnv1']['gens'][g])
            assert os.path.exists(info['cnv2']['gens'][g])

            with mdtraj.open(info['cnv1']['gens'][g]) as tfile:
                with mdtraj.open(info['cnv2']['gens'][g]) as dcdfile:
                    xyz, time, step, box = tfile.read()
                    print("Shape", xyz.shape)
                    assert xyz.shape == (11, 22, 3), xyz.shape

                    xyz_nc, time, lengths, angles = dcdfile.read()
                    np.testing.assert_array_equal(xyz * 10, xyz_nc)


@with_setup(generate_bw, cleanup)
def test_cnv_bw():
    project = Project("v1", "data/v1", 'bw')
    raw_infos = list(project.get_infos())
    nfo_infos = list(map(record(project.nfo), raw_infos))
    cnv1_infos = list(map(record(project.cnv1), nfo_infos))
    cnv2_infos = list(map(record(project.cnv2), cnv1_infos))

    assert len(cnv2_infos) > 0

    for info in cnv2_infos:
        assert os.path.exists(info['cnv2']['gens'][0])
        assert "/cnv2" in info['cnv2']['outdir']

        with mdtraj.open(info['cnv2']['gens'][0]) as tfile:
            xyz, time, step, box = tfile.read()
            print("Shape", xyz.shape)
            assert xyz.shape == (11, 22, 3), xyz.shape


@with_setup(mock2, cleanup2)
def test_cnv_nc():
    os.remove("processed/p9761/24/7/cnv.nc")
    with open("processed/p9761/24/7/info.json") as f:
        info = json.load(f)
    process._cnv2(info, chunk=2)

    top = 'tops-p9712/{top[struct]}.prmtop'.format(**info)
    trj1 = mdtraj.load(info['cnv1']['gens'][0], top=top)
    trj2 = mdtraj.load(info['cnv2']['gens'][0], top=top)

    np.testing.assert_array_almost_equal(trj1.xyz, trj2.xyz)
    np.testing.assert_array_almost_equal(trj1.unitcell_vectors,
                                         trj2.unitcell_vectors)
    np.testing.assert_array_almost_equal(trj1.time, trj2.time)


@with_setup(mock2, cleanup2)
def test_nc_cpptraj():
    os.remove("processed/p9761/24/7/cnv.nc")
    with open("processed/p9761/24/7/info.json") as f:
        info = json.load(f)
    process._cnv2(info, chunk=2)

    top = 'tops-p9712/{top[struct]}.prmtop'.format(**info)
    out = "{workdir}/test.nc".format(**info['path'])
    subprocess.check_call([
        'cpptraj',
        '-p', top,
        '-y', info['cnv2']['gens'][0],
        '-x', out
    ], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    trj1 = mdtraj.load(out, top=top)
    trj2 = mdtraj.load(info['cnv2']['gens'][0], top=top)

    np.testing.assert_array_almost_equal(trj1.xyz, trj2.xyz)
    np.testing.assert_array_almost_equal(trj1.unitcell_vectors,
                                         trj2.unitcell_vectors)
    # np.testing.assert_array_almost_equal(trj1.time, trj2.time)
