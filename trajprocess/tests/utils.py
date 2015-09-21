"""Tools for setting up a fake directory structure for processing."""

from tempfile import mkdtemp
import os
import shutil

import numpy as np
from mdtraj.formats import XTCTrajectoryFile


def write_traj(path, i):
    n_frame = 10
    with XTCTrajectoryFile(path, 'w') as f:
        f.write(np.random.randn(n_frame, 7, 3),
                time=np.arange(n_frame) + n_frame * i)


def write_run_clone(proj, run, clone, gens=None):
    if gens is None:
        gens = [0, 1]

    rc = "data/PROJ{proj}/RUN{run}/CLONE{clone}/".format(proj=proj, run=run,
                                                         clone=clone)
    os.makedirs(rc, exist_ok=True)
    for gen in gens:
        write_traj("{}/frame{}.xtc".format(rc, gen), gen)


def generate_project():
    global wd
    wd = mkdtemp()
    os.chdir(wd)
    write_run_clone(1234, 5, 7)
    write_run_clone(1234, 6, 0)


def cleanup():
    shutil.rmtree(wd)
