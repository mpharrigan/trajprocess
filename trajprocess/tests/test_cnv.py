import re
import subprocess
import os

from nose import with_setup
import mdtraj

from .utils import generate_project, cleanup, write_run_clone
from trajprocess.files import Project, record


@with_setup(generate_project, cleanup)
def test_cnv():
    project = Project("p1234", "data/PROJ1234", 'xa4')
    raw_infos = list(project.get_infos())
    nfo_infos = list(map(record(project.nfo), raw_infos))
    cat_infos = list(map(record(project.cat), nfo_infos))
    cnv_infos = list(map(record(project.cnv), cat_infos))

    assert len(cnv_infos) > 0
    return True

    for info in cat_infos:
        assert info['cat']['gen'] == 2  # per construction
        assert info['cat']['success']
        assert os.path.exists(info['cat']['xtc_out'])

        with mdtraj.open(info['cat']['xtc_out']) as tfile:
            xyz, time, step, box = tfile.read()
            print("Shape", xyz.shape)
            assert xyz.shape == (20, 7, 3)
