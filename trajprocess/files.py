"""High-level machinery for running multiple processing steps."""

from multiprocessing import Pool
import glob
import json
import os
import logging
from pathlib import Path

from . import process

log = logging.getLogger(__name__)


class Project:
    def __init__(self, code, indir, mdtype):
        assert mdtype in ['x21', 'xa4', 'bw']
        assert indir[-1] != '/', "Don't include the trailing slash!"

        self.code = code
        self.indir = indir
        self.mdtype = mdtype

        if mdtype == 'x21':
            self.nfo = process.nfo_21
            self.cat = process.cat_21
            self.cnv = process.cnv_21
        elif mdtype == 'xa4':
            self.nfo = process.nfo_a4
            self.cat = process.cat_a4
            self.cnv = process.cnv_a4
        elif mdtype == 'bw':
            pass

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

        return ([
                    {'raw': {'indir': indir,
                             'real_indir': os.path.realpath(indir)},
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
    for project in projects:
        log.info("Starting project {}".format(project))
        raw_infos = list(project.get_infos())
        log.debug("Found {} infos".format(len(raw_infos)))
        with Pool() as pool:
            nfo_infos = pool.map(record(project.nfo), raw_infos)
            cat_infos = pool.map(record(project.cat), nfo_infos, chunksize=1)
            cnv_infos = pool.map(record(project.cnv), cat_infos, chunksize=1)


def main_nav():
    return process_projects(
        Project('p9704', 'data/PROJ9704', 'x21'),
        Project('p9752', 'data/PROJ9752', 'xa4'),
        # Project('v4', 'data/v4', 'bw'),
        # Project('v5', 'data/v5', 'bw'),
    )


def main_trek():
    return process_projects(
        Project('p9712', 'PROJ9712', 'x21'),
        Project('p9761', 'PROJ9761', 'xa4'),
    )
