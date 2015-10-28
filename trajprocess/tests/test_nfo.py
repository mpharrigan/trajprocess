import os
from multiprocessing import Pool
import json

from nose import with_setup

from trajprocess.files import Processor, Trajectory, _record
from .mock1 import generate_project, cleanup, generate_bw


def _process_trajectory(trajectory):
    trajectory.info = _record(trajectory.processor.nfo, trajectory.info)
    return trajectory.info


@with_setup(generate_project, cleanup)
def test_trajectory_nfo():
    proj = Processor("p1234", "data/PROJ1234", 'xa4')
    trajectories = [Trajectory(info, proj, None) for info in proj.get_infos()]

    with Pool() as pool:
        nfo_infos = pool.map(_process_trajectory, trajectories, chunksize=1)

    n = 0
    for info in nfo_infos:
        assert os.path.exists(info['path']['info'])
        with open(info['path']['info']) as f:
            reconstitute = json.load(f)
        assert reconstitute == info
        n += 1
    assert n > 0


@with_setup(generate_project, cleanup)
def test_nfo():
    proj = Processor("p1234", "data/PROJ1234", 'xa4')
    trajectories = [Trajectory(info, proj, None) for info in proj.get_infos()]

    with Pool() as pool:
        nfo_infos = pool.map(_process_trajectory, trajectories, chunksize=1)

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
    proj = Processor("v1", "data/v1", 'bw')
    trajectories = [Trajectory(info, proj, None) for info in proj.get_infos()]

    with Pool() as pool:
        nfo_infos = pool.map(_process_trajectory, trajectories, chunksize=1)

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
