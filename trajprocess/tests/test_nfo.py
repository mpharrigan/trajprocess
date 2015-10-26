import os
from multiprocessing import Pool

from nose import with_setup

from trajprocess.files import Project, record
from .mock1 import generate_project, cleanup, generate_bw


@with_setup(generate_project, cleanup)
def test_nfo():
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
        assert set(info['raw'].keys()) == {'indir', 'real_indir', 'date',
                                           'initdate', 'gens', 'gen_glob',
                                           'success'}
        assert set(info['meta'].keys()) == {'project', 'run', 'clone'}
        assert set(info['path'].keys()) == {'info', 'workdir'}
        assert set(info['top'].keys()) == {'struct', 'fext'}


@with_setup(generate_bw, cleanup)
def test_nfo_bw():
    project = Project("v1", "data/v1", 'bw')
    raw_infos = list(project.get_infos())
    with Pool() as pool:
        nfo_infos = pool.map(record(project.nfo), raw_infos)

    assert len(nfo_infos) == 2

    for info in nfo_infos:
        assert info['meta']['project'] == 'v1'

        if info['meta']['run'] == 5:
            assert info['meta']['clone'] == 0
        elif info['meta']['run'] == 6:
            assert info['meta']['clone'] == 0
        else:
            assert False

        assert os.path.exists(info['path']['workdir'])

        assert set(info.keys()) == {'raw', 'meta', 'path', 'top'}, info
        assert set(info['raw'].keys()) == {'indir', 'real_indir', 'date',
                                           'initdate', 'gens', 'gen_glob',
                                           'success'}
        assert set(info['meta'].keys()) == {'project', 'run', 'clone'}
        assert set(info['path'].keys()) == {'info', 'workdir'}
        assert set(info['top'].keys()) == {'struct', 'fext'}
