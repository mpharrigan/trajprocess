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


def write_run_clone(wd, proj, run, clone):
    rc = "{wd}/data/PROJ{proj}/RUN{run}/CLONE{clone}/".format(wd=wd, proj=proj,
                                                              run=run,
                                                              clone=clone)
    os.makedirs(rc)
    write_traj("{}/frame0.xtc".format(rc), 0)
    write_traj("{}/frame1.xtc".format(rc), 1)


def generate_project():
    global wd
    wd = mkdtemp()
    write_run_clone(wd, 1234, 5, 7)
    write_run_clone(wd, 1234, 6, 0)
    os.chdir(wd)


def cleanup():
    shutil.rmtree(wd)
