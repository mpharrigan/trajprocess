from trajprocess import main_nav
from nose import with_setup
from .mock4 import mock_project as mock4, cleanup as cleanup4
import os
import shutil
import subprocess
import numpy as np

import json
from pathlib import Path

import mdtraj

PROJ9704_FRAMES_PER_GEN = 16
PROJ9752_FRAMES_PER_GEN = 8


# Please run ./get_addtl_ref_data.sh before attempting these tests


def _nav_asserts1(I):
    out = Path("processed.v2")
    assert out.exists()
    with (out / "p9704/9/19/info.json").open('r') as file:
        info = json.load(file)

    top = mdtraj.load_prmtop(info['stp']['outtop'])

    for step in ['raw', 'cnv2', 'stp', 'ctr']:
        assert len(info[step]['gens']) == I, (step, I)

    for i in range(I):
        assert Path(info['ctr']['gens'][i]).exists(), i

    # Load all gens
    traj = mdtraj.load(info['ctr']['gens'], top=top)
    assert len(traj) == traj.n_frames
    assert traj.n_atoms == 30962, traj.n_atoms
    assert traj.n_frames == PROJ9704_FRAMES_PER_GEN * I, traj.n_frames


def _nav_asserts2(I):
    out = Path("processed.v2")
    assert out.exists()
    with (out / "p9752/88/0/info.json").open('r') as file:
        info = json.load(file)

    assert info['cnv2']['had_overlapping_frames']

    top = mdtraj.load_prmtop(info['stp']['outtop'])

    for step in ['raw', 'cnv1', 'cnv2', 'stp', 'ctr']:
        assert len(info[step]['gens']) == I, (step, I)

    for i in range(I):
        assert Path(info['ctr']['gens'][i]).exists(), i

    # Load all gens
    traj = mdtraj.load(info['ctr']['gens'], top=top)
    assert len(traj) == traj.n_frames
    assert traj.n_atoms == 30962, traj.n_atoms
    assert traj.n_frames == PROJ9752_FRAMES_PER_GEN * I, traj.n_frames


@with_setup(mock4, cleanup4)
def test_nav():
    # Make a concatenated
    subprocess.check_call(
        ['gmx', 'trjcat', '-f'] + ['PROJ9752/RUN88/CLONE0/frame{}.xtc'.format(i)
                                   for i in range(3)] + ['-o', 'catty.xtc'],
        stderr=subprocess.STDOUT, stdout=subprocess.DEVNULL)

    # Hide one trajectory
    fromplace1 = "PROJ9704/RUN9/CLONE19/results-002"
    fromplace2 = "PROJ9752/RUN88/CLONE0/frame2.xtc"
    tmpplace1 = "hide-me1"
    tmpplace2 = "hide-me2"
    shutil.move(fromplace1, tmpplace1)
    shutil.move(fromplace2, tmpplace2)

    # Do with hidden
    main_nav()
    _nav_asserts1(2)
    _nav_asserts2(2)

    # Bring back hidden
    shutil.move(tmpplace1, fromplace1)
    assert os.path.exists(fromplace1 + "/positions.xtc")
    main_nav()
    _nav_asserts1(3)
    _nav_asserts2(2)

    # Bring back other hidden
    shutil.move(tmpplace2, fromplace2)
    assert os.path.exists(fromplace2)
    main_nav()
    _nav_asserts1(3)
    _nav_asserts2(3)
