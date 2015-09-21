import re
import subprocess
import os

from nose import with_setup
import mdtraj

from .utils import generate_project, cleanup, write_run_clone
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
    assert int(rev) >= 6
    # There is a bug in 5.0.x < 6 where appending with trjcat results in
    # segfault. http://redmine.gromacs.org/issues/1705


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
            assert xyz.shape == (20, 22, 3)


@with_setup(generate_project, cleanup)
def test_cat_from_existing():
    project = Project("p1234", "data/PROJ1234", 'xa4')
    raw_infos = list(project.get_infos())
    nfo_infos = list(map(record(project.nfo), raw_infos))
    cat_infos = list(map(record(project.cat), nfo_infos))
    assert len(cat_infos) > 0

    write_run_clone(1234, 5, 7, gens=[2])
    os.remove("data/PROJ1234/RUN5/CLONE7/frame0.xtc")
    os.remove("data/PROJ1234/RUN5/CLONE7/frame1.xtc")
    with open("data/PROJ1234/RUN5/CLONE7/frame0.xtc", 'w') as f:
        f.write("dummy!")
    with open("data/PROJ1234/RUN5/CLONE7/frame1.xtc", 'w') as f:
        f.write("dummy!")
    project = Project("p1234", "data/PROJ1234", 'xa4')
    raw_infos = list(project.get_infos())
    nfo_infos = list(map(record(project.nfo), raw_infos))
    cat_infos = list(map(record(project.cat), nfo_infos))
    assert len(cat_infos) > 0

    found_one = False
    for info in cat_infos:
        if info['meta']['run'] == 5:
            found_one = True
            assert info['cat']['gen'] == 3
            with mdtraj.open(info['cat']['xtc_out']) as tfile:
                xyz, time, step, box = tfile.read()
                print("Shape", xyz.shape)
                assert xyz.shape == (30, 22, 3), xyz.shape
        else:
            assert info['cat']['gen'] == 2
            with mdtraj.open(info['cat']['xtc_out']) as tfile:
                xyz, time, step, box = tfile.read()
                print("Shape", xyz.shape)
                assert xyz.shape == (20, 22, 3), xyz.shape
    assert found_one
