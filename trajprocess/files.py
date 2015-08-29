__author__ = 'harrigan'

from multiprocessing import Pool
import glob
import json
from . import process


class Project:
    def __init__(self, code, indir, mdtype):
        assert mdtype in ['x21, xa4, bw']

        self.code = code
        self.indir = indir
        self.mdtype = mdtype

        if mdtype == 'x21':
            self.nfo = process.nfo_21
            self.cat = process.cat_21
            self.cnv = process.cnv_21
        elif mdtype == 'a4':
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


def write_infos(infos):
    for info in infos:
        with open(info['nfo_nfoout'], 'w') as f:
            json.dump(info, f, indent=2)


def main(*projects):
    for project in projects:
        run_clone_dirs = project.get_run_clone_dirs()
        with Pool() as pool:
            nfo_infos = pool.imap_unordered(project.nfo, run_clone_dirs)
            cat_infos = pool.imap_unordered(project.cat, nfo_infos)
            cnv_infos = pool.imap_unordered(project.cnv, cat_infos)

            write_infos(cnv_infos)


def sample_main():
    return main(
        Project('p9704', 'data/PROJ9704', 'x21'),
        Project('p9751', 'data/PROJ9751', 'xa4'),
        # Project('v4', 'data/v4', 'bw'),
        # Project('v5', 'data/v5', 'bw'),
    )
