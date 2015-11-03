"""Functions for performing the individual processing steps.

 - "nfo": Prepare meta-information for each trajectory
 - "cnv1": Run trjconv for pbc imaging (gromacs runs only)
 - "cnv2": Convert files to netcdf

"""

import subprocess
import os
import glob
import re
import json
import logging
from operator import itemgetter
from datetime import datetime

from mdtraj.formats import XTCTrajectoryFile, NetCDFTrajectoryFile
import mdtraj.utils
import numpy as np

log = logging.getLogger(__name__)


class config:
    prefix = "processed.v2"


def _nfo(info, *, rncln_re, gen_glob, gen_re, gen=None, clone=None):
    if 'run' not in info['meta']:
        # Get metadata
        rncln_ma = rncln_re.search(info['raw']['indir'])
        info['meta']['run'] = int(rncln_ma.group(1))
        info['meta']['clone'] = (int(rncln_ma.group(2))
                                 if clone is None else clone)
        log.debug("Got metadata {meta[project]}-{meta[run]}-{meta[clone]}"
                  .format(**info))

    if 'path' not in info:
        path = {'workdir': "{prefix}/{project}/{run}/{clone}"
            .format(prefix=config.prefix, **info['meta'])}
        path['info'] = "{workdir}/info.json".format(**path)
        info['path'] = path
        os.makedirs(info['path']['workdir'], exist_ok=True)
        log.debug("Make workdir: {path[workdir]}".format(**info))

    # Get gens
    raw = info['raw']
    raw['gen_glob'] = gen_glob
    raw['date'] = datetime.now().isoformat()
    raw['gens'] = []  # re-do each time

    gens = sorted(
        ((int(gen_re.search(gen_fn).group(1) if gen is None else gen), gen_fn)
         for gen_fn in glob.glob("{indir}/{gen_glob}".format(**raw))),
        key=itemgetter(0)
    )

    # Make sure they're contiguous
    prev_gen = -1
    for gen, gen_fn in gens:
        raw['gens'] += [gen_fn]

        if gen != prev_gen + 1:
            log.error("Found discontinous gens "
                      "in {meta[project]}-{meta[run]}-{meta[clone]}. "
                      "It went from {i1} to {i2}."
                      .format(i1=prev_gen, i2=gen, **info))
            raw['success'] = False
            info['raw'] = raw
            return info
        prev_gen = gen

    if 'exclude' in raw:
        log.warning("Excluding {meta[project]}-{meta[run]}-{meta[clone]}. "
                    "Reason: {raw[exclude]}".format(**info))
        raw['gens'] = []
        raw['success'] = False
    else:
        raw['success'] = True

    info['raw'] = raw

    # Get structure (topology) data
    if 'top' not in info:
        struct_fn = "structs-{meta[project]}.json".format(**info)
        try:
            with open(struct_fn) as f:
                stru = json.load(f)
                string_key = str(info['meta']['run'])  # ugh
                info['top'] = stru[string_key]
        except Exception as e:
            log.warning("No structure information. {}".format(e))

    log.info("NFO: {project} run {run} clone {clone}".format(**info['meta']))

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


def _run_trjconv(info, gen, gen_fn):
    out_fn = "{outdir}/{gen}.{outext}".format(gen=gen, **info['cnv1'])
    log.debug("Running trjconv {} {}".format(gen_fn, out_fn))
    with open(info['cnv1']['log'], 'a') as logf:
        popen = subprocess.Popen([
            'gmx', 'trjconv',
            '-f', gen_fn,
            '-o', out_fn,
            '-s', info['cnv1']['topology'],
            '-pbc', 'mol',
            '-center',
            '-skip', str(info['cnv1']['stride']),
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
    return out_fn


def _cnv1(info, *, stride, topology, skip=False):
    info['cnv1'] = {
        'date': datetime.now().isoformat(),
        'stride': stride,
        'topology': topology,
        'skip': skip,
        'log': "{workdir}/cnv1.log".format(**info['path']),
        'outdir': "{workdir}/cnv1".format(**info['path']),
        'outext': 'xtc',
        'gens': [] if 'cnv1' not in info else info['cnv1']['gens'],
    }

    if not info['raw']['success']:
        info['cnv1']['success'] = False
        return info

    if skip:
        info['cnv1']['success'] = True
        return info

    os.makedirs(info['cnv1']['outdir'], exist_ok=True)

    done = len(info['cnv1']['gens'])
    log.info("CNV1: {meta[project]}-{meta[run]}-{meta[clone]}. "
             "Done {done}, doing {todo}"
             .format(done=done, todo=len(info['raw']['gens']) - done, **info))
    for gen, gen_fn in enumerate(info['raw']['gens']):
        if gen < done:
            continue
        out_fn = _run_trjconv(info, gen, gen_fn)
        info['cnv1']['gens'] += [out_fn]

    info['cnv1']['success'] = True
    return info


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


def _nc_a_traj(info, gen, gen_fn, has_overlapping_frames):
    out_fn = "{outdir}/{gen}.{outext}".format(gen=gen, **info['cnv2'])
    log.debug("Converting to netcdf {} {}".format(gen_fn, out_fn))
    with XTCTrajectoryFile(gen_fn, 'r') as xtc:
        with NetCDFTrajectoryFile(out_fn, 'w') as nc:
            _nc_a_chunk(xtc, nc, has_overlapping_frames)

    return out_fn


def _cnv2(info, *, has_overlapping_frames, chunk=100):
    info['cnv2'] = {
        'date': datetime.now().isoformat(),
        'chunk': chunk,
        'had_overlapping_frames': has_overlapping_frames,
        'log': "{workdir}/cnv2.log".format(**info['path']),
        'outdir': "{workdir}/cnv2".format(**info['path']),
        'outext': 'nc',
        'gens': [] if 'cnv2' not in info else info['cnv2']['gens'],
    }

    if not info['cnv1']['success']:
        info['cnv2']['success'] = False
        return info

    if info['cnv1']['skip']:
        prev_gens = info['raw']['gens']
    else:
        prev_gens = info['cnv1']['gens']

    done = len(info['cnv2']['gens'])
    log.info("CNV2: {meta[project]}-{meta[run]}-{meta[clone]}. "
             "Converting to nc. Done {done}, todo {todo}"
             .format(done=done, todo=len(prev_gens) - done, **info))
    os.makedirs(info['cnv2']['outdir'], exist_ok=True)
    for gen, gen_fn in enumerate(prev_gens):
        if gen < done:
            continue
        out_fn = _nc_a_traj(info, gen, gen_fn, has_overlapping_frames)
        info['cnv2']['gens'] += [out_fn]

    info['cnv2']['success'] = True
    return info


def cnv1(info, projcode):
    if info['meta']['project'] == 'p9752':
        stride = 4
    elif info['meta']['project'] == 'p9761':
        stride = 8
    else:
        stride = 1

    topology = None
    skip = False
    if projcode == 'xa4':
        topology = "{raw[indir]}/frame0.tpr".format(**info)
    elif projcode == 'bw':
        topology = "{raw[indir]}/topol.tpr".format(**info)
    else:
        skip = True

    return _cnv1(
        info,
        stride=stride,
        topology=topology,
        skip=skip,
    )


def cnv2(info, projcode):
    if projcode == 'xa4':
        overlap = True
    else:
        overlap = False

    return _cnv2(
        info,
        has_overlapping_frames=overlap
    )
