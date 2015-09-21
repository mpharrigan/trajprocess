"""Tools for setting up a fake directory structure for processing."""

from tempfile import mkdtemp
import os
import shutil
import json

from pkg_resources import resource_filename


# command for generating reference data:
# gmx mdrun -nsteps 5000 -s frame0.tpr -cpi -noappend
#
# Do that three times.


def write_run_clone(proj, run, clone, gens=None):
    if gens is None:
        gens = [0, 1]

    rc = "data/PROJ{proj}/RUN{run}/CLONE{clone}/".format(proj=proj, run=run,
                                                         clone=clone)
    os.makedirs(rc, exist_ok=True)
    tpr_fn = resource_filename(__name__, 'topol.tpr')
    shutil.copy(tpr_fn, "{}/frame0.tpr".format(rc))
    for gen in gens:
        shutil.copy(resource_filename(__name__,
                                      "traj_comp.part{:04d}.xtc".format(
                                          gen + 1)),
                    "{}/frame{}.xtc".format(rc, gen))


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
