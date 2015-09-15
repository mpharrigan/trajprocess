__author__ = 'harrigan'

from tempfile import mkdtemp

import numpy as np
import os
from mdtraj.formats import XTCTrajectoryFile
import shutil


def write_traj(path):
    with XTCTrajectoryFile(path, 'w') as f:
        f.write(np.random.randn(10, 7, 3))


def write_run_clone(wd, proj, run, clone):
    rc = "{wd}/data/PROJ{proj}/RUN{run}/CLONE{clone}/".format(wd=wd, proj=proj,
                                                              run=run,
                                                              clone=clone)
    os.makedirs(rc)
    write_traj("{}/frame0.xtc".format(rc))
    write_traj("{}/frame1.xtc".format(rc))


def generate_project():
    global wd
    wd = mkdtemp()
    write_run_clone(wd, 1234, 5, 7)
    write_run_clone(wd, 1234, 6, 0)
    os.chdir(wd)


def cleanup():
    shutil.rmtree(wd)
