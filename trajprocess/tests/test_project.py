import os
from multiprocessing import Pool
import json
import logging

from nose import with_setup

from trajprocess.files import Processor, _record
from .mock1 import generate_project, cleanup

logging.basicConfig(level=logging.DEBUG)


@with_setup(generate_project, cleanup)
def test_glob():
    proj = Processor("p1234", "data/PROJ1234", 'xa4')
    assert len(proj.get_run_clone_dirs()) > 0


@with_setup(generate_project, cleanup)
def test_infogen():
    proj = Processor("p1234", "data/PROJ1234", 'xa4')
    infos = list(proj.get_infos())
    assert len(infos) > 0, "{}".format(infos)


@with_setup(generate_project, cleanup)
def test_infogen_2():
    proj = Processor("p1234", "data/PROJ1234", 'xa4')
    n = 0
    for info in proj.get_infos():
        assert os.path.exists(info['raw']['indir'])
        n += 1

    assert n > 0


@with_setup(generate_project, cleanup)
def test_nfo_1():
    proj = Processor("p1234", "data/PROJ1234", 'xa4')
    raw_infos = proj.get_infos()
    nfo_infos = [proj.nfo(ri) for ri in raw_infos]
    assert len(nfo_infos) > 0
    for info in nfo_infos:
        assert info['raw']['success']


@with_setup(generate_project, cleanup)
def test_nfo_record():
    proj = Processor("p1234", "data/PROJ1234", 'xa4')
    raw_infos = proj.get_infos()
    nfo_infos = [_record(proj.nfo, ri) for ri in raw_infos]
    assert len(nfo_infos) > 0

    for info in nfo_infos:
        assert os.path.exists(info['path']['info'])
        with open(info['path']['info']) as f:
            reconstitute = json.load(f)
        assert reconstitute == info


@with_setup(generate_project, cleanup)
def test_nfo_pool():
    proj = Processor("p1234", "data/PROJ1234", 'xa4')
    raw_infos = proj.get_infos()
    print("Raw infos", list(proj.get_infos()))
    with Pool() as pool:
        print("Function", proj.nfo)
        nfo_infos = pool.map(proj.nfo, raw_infos)

    print("Returned", nfo_infos)
    assert len(list(nfo_infos)) > 0


@with_setup(generate_project, cleanup)
def test_tpr():
    assert os.path.exists("data/PROJ1234/RUN5/CLONE7/frame0.tpr")
