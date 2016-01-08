"""Functions for performing the individual processing steps.

 - "nfo": Prepare meta-information for each trajectory
 - "cnv1": Run trjconv for pbc imaging (gromacs runs only)
 - "cnv2": Convert files to netcdf

"""

import logging
import re
import subprocess
from datetime import datetime

import mdtraj.utils
import numpy as np
from mdtraj.formats import XTCTrajectoryFile, NetCDFTrajectoryFile

log = logging.getLogger(__name__)


class config:
    prefix = "processed.v2"


def _nfo(info, *, rncln_re, gen_glob, gen_re, gen=None, clone=None):
    return info


def nfo(info, projcode):
    rncln_res = {
        'xa4': re.compile(r"RUN(\d+)/CLONE(\d+)/"),
        'x21': re.compile(r"RUN(\d+)/CLONE(\d+)/"),
        'bw': re.compile(r"run-(\d+)/"),
    }
    gen_globs = {
        'xa4': "frame*.xtc",
        'x21': "results-???/positions.xtc",
        'bw': "traj_comp.xtc",
    }
    gen_res = {
        'xa4': re.compile(r"frame(\d+).xtc"),
        'x21': re.compile(r"results-(\d+)/positions.xtc"),
        'bw': re.compile(r""),
    }
    gen = None
    clone = None
    if projcode == 'bw':
        gen = 0
        clone = 0

    return _nfo(
            info,
            rncln_re=rncln_res[projcode],
            gen_glob=gen_globs[projcode],
            gen_re=gen_res[projcode],
            gen=gen,
            clone=clone,
    )


def run_trjconv(in_fn, out_fn, *, log_fn, top_fn, stride):
    with open(log_fn, 'a') as logf:
        popen = subprocess.Popen([
            'gmx', 'trjconv',
            '-f', in_fn,
            '-o', out_fn,
            '-s', top_fn,
            '-pbc', 'mol',
            '-center',
            '-skip', str(stride),
        ],
                stdin=subprocess.PIPE,
                stdout=logf,
                stderr=subprocess.STDOUT
        )
        # Center based on 1 - Protein
        # Output 0 - System
        popen.communicate(b"1\n0")
        popen.wait()

        if popen.returncode != 0:
            raise RuntimeError("Non-zero exit code from trjconv {}"
                               .format(popen.returncode))


def _nc_a_chunk(xtc, nc, has_overlapping_frames):
    xyz, time, step, box = xtc.read()
    assert box.ndim == 3, box.ndim
    al, bl, cl, alpha, beta, gamma = \
        mdtraj.utils.box_vectors_to_lengths_and_angles(
                box[:, 0, :], box[:, 1, :], box[:, 2, :]
        )

    xyz = xyz * 10
    blengs = np.asarray([al, bl, cl]).T * 10
    bangles = np.asarray([alpha, beta, gamma]).T

    sl = slice(0, -1 if has_overlapping_frames else None, 1)
    nc.write(
            xyz[sl, ...],
            time[sl, ...],
            blengs[sl, ...],
            bangles[sl, ...],
    )


def convert_to_nc(in_fn, out_fn, *, has_overlapping_frames):
    with XTCTrajectoryFile(in_fn, 'r') as xtc:
        with NetCDFTrajectoryFile(out_fn, 'w') as nc:
            _nc_a_chunk(xtc, nc, has_overlapping_frames)



# TODO: remove below

def cnv1(info, projcode):
    if info['meta']['project'] == 'p9752':
        stride = 4
    elif info['meta']['project'] == 'p9761':
        stride = 8
    else:
        stride = 1

    if projcode == 'xa4':
        topology = "{raw[indir]}/frame0.tpr".format(**info)
    elif projcode == 'bw':
        topology = "{raw[indir]}/topol.tpr".format(**info)
    else:
        raise Exception

    run_trjconv(
            info['infn'], info['outfn'],
            log_fn=info['log_fn'],
            top_fn=topology,
            stride=stride,
    )


def cnv2(info, projcode):
    if projcode == 'xa4':
        overlap = True
    else:
        overlap = False

    return convert_to_nc(
            info['infn'], info['outfn'],
            has_overlapping_frames=overlap
    )
