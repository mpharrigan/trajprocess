"""High-level machinery for running multiple processing steps."""

from multiprocessing import Pool
import glob
import json
import os
import logging
from pathlib import Path
from datetime import datetime

from . import process, postprocess

log = logging.getLogger(__name__)


class _processWrap:
    def __init__(self, func, mdtype):
        self.func = func
        self.mdtype = mdtype

    def __call__(self, info):
        return self.func(info, self.mdtype)


class Project:
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

    def get_infos(self):
        prev = Path("processed/{project}".format(project=self.code))
        prev_indirs = set()
        prev_infos = list()
        for info_fn in prev.glob("**/info.json"):
            with info_fn.open() as f:
                info = json.load(f)
                prev_infos.append(info)
                prev_indirs.add(info['raw']['indir'])

        new_indirs = set(self.get_run_clone_dirs()) - prev_indirs
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
        return "{code} ({indir})".format(**self.__dict__)


class record:
    def __init__(self, func):
        self.func = func

    def __call__(self, info):
        info = self.func(info)
        with open(info['path']['info'], 'w') as f:
            json.dump(info, f, indent=2)
        return info


def process_projects(*projects):
    infos = []
    for project in projects:
        log.info("Starting project {}".format(project))
        raw_infos = list(project.get_infos())
        log.info("Found {} infos".format(len(raw_infos)))
        with Pool() as pool:
            nfo_infos = pool.map(record(project.nfo), raw_infos)
            cnv1_infos = pool.map(record(project.cnv1), nfo_infos, chunksize=1)
            cnv2_infos = pool.map(record(project.cnv2), cnv1_infos, chunksize=1)
        infos += cnv2_infos
    return infos


class _postprocessWrap:
    def __init__(self, func, systemcode):
        self.func = func
        self.systemcode = systemcode

    def __call__(self, info):
        return self.func(info, self.systemcode)


class Postprocess:
    def __init__(self, system):
        self.system = system

        self.stp = _postprocessWrap(postprocess.stp, system)
        self.ctr = _postprocessWrap(postprocess.ctr, system)


def process_post(postprocessor, cnv_infos):
    with Pool() as pool:
        stp_infos = pool.map(record(postprocessor.stp), cnv_infos, chunksize=1)
        ctr_infos = pool.map(record(postprocessor.ctr), stp_infos, chunksize=1)
    return ctr_infos


def main_nav():
    infos = process_projects(
        Project('p9704', 'PROJ9704', 'x21'),
        Project('p9752', 'PROJ9752', 'xa4'),
        Project('v4', 'v4', 'bw'),
        Project('v5', 'v5', 'bw'),
    )
    return process_post(Postprocess('nav'), infos)


def main_trek():
    infos = process_projects(
        Project('p9712', 'PROJ9712', 'x21'),
        Project('p9761', 'PROJ9761', 'xa4'),
    )
    return process_post(Postprocess('trek'), infos)
