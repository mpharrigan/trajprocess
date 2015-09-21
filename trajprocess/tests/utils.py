"""Tools for setting up a fake directory structure for processing."""

from tempfile import mkdtemp
import os
import shutil
import json

import numpy as np
from mdtraj.formats import XTCTrajectoryFile

from pkg_resources import resource_filename


def write_traj(path, i):
    n_frame = 10
    with XTCTrajectoryFile(path, 'w') as f:
        f.write(np.random.randn(n_frame, 22, 3),
                time=np.arange(n_frame) + n_frame * i)


def write_run_clone(proj, run, clone, gens=None):
    if gens is None:
        gens = [0, 1]

    rc = "data/PROJ{proj}/RUN{run}/CLONE{clone}/".format(proj=proj, run=run,
                                                         clone=clone)
    os.makedirs(rc, exist_ok=True)
    for gen in gens:
        write_traj("{}/frame{}.xtc".format(rc, gen), gen)

    tpr_fn = resource_filename(__name__, 'topol.tpr')
    shutil.copy(tpr_fn, "{}/frame0.tpr".format(rc))


def generate_project():
    global wd
    wd = mkdtemp()
    os.chdir(wd)
    write_run_clone(1234, 5, 7)
    write_run_clone(1234, 6, 0)
    with open('structs-p1234.json', 'w') as f:
        json.dump({
            5: {'struct': 'stru1', 'fext': 'pdb'},
            6: {'struct': 'stru2', 'fext': 'pdb'}
        }, f)


def cleanup():
    shutil.rmtree(wd)
