from trajprocess import main_trek
from nose import with_setup
from .mock3 import mock_project as mock3, cleanup as cleanup3
import os
import shutil
import subprocess
import numpy as np

import json
from pathlib import Path

import mdtraj

PROJ61_LENGTH_PER_GEN = 2


# Requires:
# http://web.stanford.edu/~harrigan/mock3-reference.tar.bz2


def _trek_asserts1(I):
    out = Path("processed.v2")
    assert out.exists()
    with (out / "p9712/5/32/info.json").open('r') as file:
        info = json.load(file)

    top = mdtraj.load_prmtop(info['stp']['outtop'])

    for step in ['raw', 'cnv2', 'stp', 'ctr']:
        assert len(info[step]['gens']) == I, (step, I)

    for i in range(I):
        assert Path(info['ctr']['gens'][i]).exists(), i

    # Load all gens
    traj = mdtraj.load(info['ctr']['gens'], top=top)
    assert len(traj) == traj.n_frames
    assert traj.n_atoms == 30962, (traj.n_atoms)
    assert traj.n_frames == 16 * I, traj.n_frames


def _trek_asserts2(I):
    out = Path("processed.v2")
    assert out.exists()
    with (out / "p9761/3/9/info.json").open('r') as file:
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
    assert traj.n_frames == PROJ61_LENGTH_PER_GEN * I, traj.n_frames


@with_setup(mock3, None)
def test_trek():
    # Make a concatenated
    subprocess.check_call(
        ['gmx', 'trjcat', '-f'] + ['PROJ9761/RUN3/CLONE9/frame{}.xtc'.format(i)
                                   for i in range(3)] + ['-o', 'catty.xtc'],
        stderr=subprocess.STDOUT, stdout=subprocess.DEVNULL)

    # Hide one trajectory
    fromplace1 = "PROJ9712/RUN5/CLONE32/results-002"
    fromplace2 = "PROJ9761/RUN3/CLONE9/frame2.xtc"
    tmpplace1 = "hide-me1"
    tmpplace2 = "hide-me2"
    shutil.move(fromplace1, tmpplace1)
    shutil.move(fromplace2, tmpplace2)

    # Do with hidden
    main_trek()
    _trek_asserts1(2)
    _trek_asserts2(2)

    # Bring back hidden
    shutil.move(tmpplace1, fromplace1)
    assert os.path.exists(fromplace1 + "/positions.xtc")
    main_trek()
    _trek_asserts1(3)
    _trek_asserts2(2)

    # Bring back other hidden
    shutil.move(tmpplace2, fromplace2)
    assert os.path.exists(fromplace2)
    main_trek()
    _trek_asserts1(3)
    _trek_asserts2(3)


@with_setup(mock3, cleanup3)
def test_lengths():
    num = 3
    inptrajs = ['PROJ9761/RUN3/CLONE9/frame{}.xtc'.format(i)
                for i in range(num)]
    stride = 8
    subprocess.check_call(
        ['gmx', 'trjcat', '-f'] + inptrajs + ['-o', 'catty.xtc'],
        stderr=subprocess.STDOUT, stdout=subprocess.DEVNULL)

    with mdtraj.open("catty.xtc") as xtc:
        stridelen = len(xtc) // stride
        remain = len(xtc) % stride
        assert stridelen == num * PROJ61_LENGTH_PER_GEN, (stridelen, remain)

    top = mdtraj.load_prmtop("tops-p9712/4bw5.prmtop")
    traj1 = mdtraj.load("catty.xtc", top=top)[::stride]
    # blarg! the last frame is duplicatey
    traj2 = mdtraj.load(inptrajs[0], top=top)[::stride][:-1]
    traj2 += mdtraj.load(inptrajs[1], top=top)[::stride][:-1]
    traj2 += mdtraj.load(inptrajs[2], top=top)[::stride]
    traj3 = mdtraj.load(inptrajs, top=top,
                        discard_overlapping_frames=True)[::stride]

    np.testing.assert_array_equal(traj1.xyz, traj3.xyz)
    np.testing.assert_array_equal(traj1.xyz, traj2.xyz)
