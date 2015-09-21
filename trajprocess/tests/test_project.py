import os
from multiprocessing import Pool
import json

from nose import with_setup

from trajprocess.files import Project, record
from .utils import generate_project, cleanup


@with_setup(generate_project, cleanup)
def test_glob():
    proj = Project("p1234", "data/PROJ1234", 'xa4')
    assert len(proj.get_run_clone_dirs()) > 0


@with_setup(generate_project, cleanup)
def test_infogen():
    proj = Project("p1234", "data/PROJ1234", 'xa4')
    infos = list(proj.get_infos())
    assert len(infos) > 0, "{}".format(infos)


@with_setup(generate_project, cleanup)
def test_infogen_2():
    proj = Project("p1234", "data/PROJ1234", 'xa4')
    n = 0
    for info in proj.get_infos():
        assert os.path.exists(info['raw']['indir'])
        n += 1

    assert n > 0


@with_setup(generate_project, cleanup)
def test_nfo():
    proj = Project("p1234", "data/PROJ1234", 'xa4')
    raw_infos = proj.get_infos()
    nfo_infos = [proj.nfo(ri) for ri in raw_infos]
    assert len(nfo_infos) > 0


@with_setup(generate_project, cleanup)
def test_nfo_record():
    proj = Project("p1234", "data/PROJ1234", 'xa4')
    raw_infos = proj.get_infos()
    func = record(proj.nfo)
    nfo_infos = [func(ri) for ri in raw_infos]
    assert len(nfo_infos) > 0

    for info in nfo_infos:
        assert os.path.exists(info['path']['info'])
        with open(info['path']['info']) as f:
            reconstitute = json.load(f)
        assert reconstitute == info


@with_setup(generate_project, cleanup)
def test_nfo_pool():
    proj = Project("p1234", "data/PROJ1234", 'xa4')
    raw_infos = proj.get_infos()
    print("Raw infos", list(proj.get_infos()))
    with Pool() as pool:
        print("Function", proj.nfo)
        nfo_infos = pool.map(proj.nfo, raw_infos)

    print("Returned", nfo_infos)
    assert len(list(nfo_infos)) > 0


@with_setup(generate_project, cleanup)
def test_nfo_record_pool():
    proj = Project("p1234", "data/PROJ1234", 'xa4')
    raw_infos = proj.get_infos()
    with Pool() as pool:
        nfo_infos = pool.map(record(proj.nfo), raw_infos)

    n = 0
    for info in nfo_infos:
        assert os.path.exists(info['path']['info'])
        with open(info['path']['info']) as f:
            reconstitute = json.load(f)
        assert reconstitute == info
        n += 1
    assert n > 0


@with_setup(generate_project, cleanup)
def test_nfo_step():
    project = Project("p1234", "data/PROJ1234", 'xa4')
    raw_infos = list(project.get_infos())
    with Pool() as pool:
        nfo_infos = pool.map(record(project.nfo), raw_infos)

    for info in nfo_infos:
        assert info['meta']['project'] == 'p1234'

        if info['meta']['run'] == 5:
            assert info['meta']['clone'] == 7
        elif info['meta']['run'] == 6:
            assert info['meta']['clone'] == 0
        else:
            assert False

        assert os.path.exists(info['path']['workdir'])

        assert set(info.keys()) == {'raw', 'meta', 'path', 'top'}, info
        assert set(info['raw'].keys()) == {'indir', 'real_indir'}
        assert set(info['meta'].keys()) == {'project', 'run', 'clone'}
        assert set(info['path'].keys()) == {'info', 'workdir'}
        assert set(info['top'].keys()) == {'struct', 'fext'}


@with_setup(generate_project, cleanup)
def test_tpr():
    assert os.path.exists("data/PROJ1234/RUN5/CLONE7/frame0.tpr")
