import os

from nose import with_setup
import mdtraj

from .utils import generate_project, cleanup, generate_bw
from trajprocess.files import Project, record


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

        with mdtraj.open(info['cnv']['xtc_out']) as tfile:
            xyz, time, step, box = tfile.read()
            print("Shape", xyz.shape)
            assert xyz.shape == (21, 22, 3), xyz.shape


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
