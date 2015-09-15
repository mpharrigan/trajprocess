from nose import with_setup
from .utils import generate_project, cleanup

import re
import subprocess
import os
import mdtraj

from trajprocess.files import Project, record


def test_gmx_version():
    gmx_vstring = subprocess.check_output(
        ['gmx', '-version'],
        universal_newlines=True).splitlines()[0]
    print(gmx_vstring)
    maj, min, rev = re.search(r"VERSION (\d+\.\d+\.\d+)", gmx_vstring) \
        .group(1).split(".")
    assert int(maj) == 5
    assert int(min) >= 0
    assert int(rev) >= 4


@with_setup(generate_project, cleanup)
def test_cat():
    project = Project("p1234", "data/PROJ1234", 'xa4')
    raw_infos = list(project.get_infos())
    nfo_infos = list(map(record(project.nfo), raw_infos))
    cat_infos = list(map(record(project.cat), nfo_infos))

    assert len(cat_infos) > 0

    for info in cat_infos:
        assert info['cat']['gen'] == 2  # per construction
        assert info['cat']['success']
        assert os.path.exists(info['cat']['xtc_out'])

        with mdtraj.open(info['cat']['xtc_out']) as tfile:
            xyz, time, step, box = tfile.read()
            print("Shape", xyz.shape)
            assert xyz.shape == (20, 7, 3)
