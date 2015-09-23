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


def write_run_clone_bw(proj, run, clone, gens=None):
    write_run_clone(
        gens,
        rc="data/{proj}/run-{run}/".format(proj=proj, run=run),
        tfn="{rc}/traj_comp.xtc",
    )


def write_run_clone_a4(proj, run, clone, gens=None):
    write_run_clone(
        gens,
        rc="data/PROJ{proj}/RUN{run}/CLONE{clone}/".format(proj=proj, run=run,
                                                           clone=clone),
        tfn="{rc}/frame{gen}.xtc"
    )


def write_run_clone(gens, rc, tfn):
    if gens is None:
        gens = [0, 1]

    os.makedirs(rc, exist_ok=True)
    tpr_fn = resource_filename(__name__, 'topol.tpr')
    shutil.copy(tpr_fn, "{}/frame0.tpr".format(rc))
    for gen in gens:
        shutil.copy(
            resource_filename(__name__,
                              "traj_comp.part{:04d}.xtc".format(gen + 1)),
            tfn.format(rc=rc, gen=gen)
        )


def generate_project():
    global wd
    wd = mkdtemp()
    os.chdir(wd)
    write_run_clone_a4(1234, 5, 7)
    write_run_clone_a4(1234, 6, 0)
    with open('structs-p1234.json', 'w') as f:
        json.dump({
            5: {'struct': 'stru1', 'fext': 'pdb'},
            6: {'struct': 'stru2', 'fext': 'pdb'}
        }, f)


def generate_bw():
    global wd
    wd = mkdtemp()
    os.chdir(wd)
    write_run_clone_bw('v1', 5, 7, gens=[0])
    write_run_clone_bw('v1', 6, 0, gens=[0])
    with open('structs-v1.json', 'w') as f:
        json.dump({
            5: {'struct': 'stru1', 'fext': 'pdb'},
            6: {'struct': 'stru2', 'fext': 'pdb'}
        }, f)


def cleanup():
    shutil.rmtree(wd)
