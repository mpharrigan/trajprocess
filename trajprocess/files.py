"""High-level machinery for running multiple processing steps."""

from multiprocessing import Pool
import glob
import json
import os
import logging
from pathlib import Path
from datetime import datetime

from . import process, postprocess
from .process import config

log = logging.getLogger(__name__)


class _processWrap:
    def __init__(self, func, which):
        self.func = func
        self.which = which

    def __call__(self, info):
        return self.func(info, self.which)


class Processor:
    def __init__(self, code, indir, mdtype):
        assert mdtype in ['x21', 'xa4', 'bw']
        assert indir[-1] != '/', "Don't include the trailing slash!"

        self.code = code
        self.indir = indir
        self.mdtype = mdtype

        self.nfo = _processWrap(process.nfo, mdtype)
        self.cnv1 = _processWrap(process.cnv1, mdtype)
        self.cnv2 = _processWrap(process.cnv2, mdtype)

    def get_run_clone_dirs(self):
        if self.mdtype in ['x21', 'xa4']:
            return glob.glob("{}/RUN*/CLONE*/".format(self.indir))
        elif self.mdtype == 'bw':
            return glob.glob("{}/run-*/".format(self.indir))
        else:
            raise ValueError

    def get_infos(self):
        prev = Path("{prefix}/{project}"
                    .format(prefix=config.prefix, project=self.code))
        prev_indirs = set()
        prev_infos = list()
        for info_fn in prev.glob("**/info.json"):
            with info_fn.open() as f:
                info = json.load(f)
                prev_infos.append(info)
                prev_indirs.add(info['raw']['indir'])

        new_indirs = set(self.get_run_clone_dirs()) - prev_indirs
        log.info("Getting infos for {}".format(self))
        log.info("Found {} previous directories".format(len(prev_indirs)))
        log.info("Found {} new directories".format(len(new_indirs)))

        return ([{'raw': {'indir': indir,
                          'real_indir': os.path.realpath(indir),
                          'initdate': datetime.now().isoformat()
                          },
                  'meta': {'project': self.code},
                  }
                 for indir in new_indirs
                 ]
                + prev_infos)

    def __repr__(self):
        return "Project {code} ({indir})".format(**self.__dict__)


class Postprocessor:
    def __init__(self, system):
        self.system = system

        self.stp = _processWrap(postprocess.stp, system)
        self.ctr = _processWrap(postprocess.ctr, system)


class Trajectory:
    def __init__(self, info, processor, postprocessor):
        self.info = info
        self.processor = processor
        self.postprocessor = postprocessor


def _record(func, info):
    info = func(info)
    with open(info['path']['info'], 'w') as f:
        json.dump(info, f, indent=2)
    return info


def _process_trajectory(trajectory):
    trajectory.info = _record(trajectory.processor.nfo, trajectory.info)
    trajectory.info = _record(trajectory.processor.cnv1, trajectory.info)
    trajectory.info = _record(trajectory.processor.cnv2, trajectory.info)
    trajectory.info = _record(trajectory.postprocessor.stp, trajectory.info)
    trajectory.info = _record(trajectory.postprocessor.ctr, trajectory.info)


def process_trajectories(*processors, postprocessor):
    trajectories = []
    for proc in processors:
        for info in proc.get_infos():
            trajectories += [Trajectory(info, proc, postprocessor)]

    with Pool() as pool:
        pool.map(_process_trajectory, trajectories, chunksize=1)

    log.info("Done!")


def main_nav():
    process_trajectories(
        Processor('p9704', 'PROJ9704', 'x21'),
        Processor('p9752', 'PROJ9752', 'xa4'),
        Processor('v4', 'v4', 'bw'),
        Processor('v5', 'v5', 'bw'),
        postprocessor=Postprocessor('nav')
    )


def main_trek():
    process_trajectories(
        Processor('p9712', 'PROJ9712', 'x21'),
        Processor('p9761', 'PROJ9761', 'xa4'),
        postprocessor=Postprocessor('trek')
    )
