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


def _nfo(info, *, rncln_re, gen_glob, gen_re, gen=None, clone=None):
    if 'meta' not in info:
        # Get metadata
        rncln_ma = rncln_re.search(info['raw']['indir'])
        meta = {
            'project': info['meta']['project'],
            'run': int(rncln_ma.group(1)),
            'clone': int(rncln_ma.group(2)) if clone is None else clone,
        }
        info['meta'] = meta

    if 'path' not in info:
        path = {'workdir': "processed/{project}/{run}/{clone}"
            .format(**info['meta'])}
        path['info'] = "{workdir}/info.json".format(**path)
        info['path'] = path

    # Get gens
    raw = info['raw']
    raw['gen_glob'] = gen_glob
    raw['date'] = datetime.now().isoformat()
    if 'gens' not in raw:
        raw['gens'] = []

    gens = sorted(
        ((int(gen_re.search(gen_fn).group(1) if gen is None else gen), gen_fn)
         for gen_fn in glob.glob("{indir}/{gen_glob}".format(**raw))),
        key=itemgetter(0)
    )

    # Make sure they're contiguous
    prev_gen = -1
    for gen, gen_fn in gens:
        raw['gens'] += [{'gen': gen, 'gen_fn': gen_fn}]

        if gen != prev_gen + 1:
            raw['success'] = False
            info['raw'] = raw
            return info

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

    # Set up working directory
    os.makedirs(info['path']['workdir'], exist_ok=True)
    log.debug("NFO: {project} run {run} clone {clone}".format(**info['meta']))

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

    log.debug("CNV1: {meta[project]}-{meta[run]}-{meta[clone]}. "
              "Starting conversion with trjconv"
              .format(**info))

    os.makedirs(info['cnv1']['outdir'], exist_ok=True)

    done = set(gen for gen, _ in info['cnv1']['gens'])
    assert len(done) == len(info['cnv1']['gens'])
    for gen, gen_fn in info['raw']['gens']:
        if gen in done:
            continue
        out_fn = _run_trjconv(info, gen, gen_fn)
        info['cnv1']['gens'] += [(gen, out_fn)]

    info['cnv1']['success'] = True
    return info


def _nc_a_chunk(xtc, nc, chunk):
    xyz, time, step, box = xtc.read(chunk)
    assert box.ndim == 3, box.ndim
    al, bl, cl, alpha, beta, gamma = \
        mdtraj.utils.box_vectors_to_lengths_and_angles(
            box[:, 0, :], box[:, 1, :], box[:, 2, :]
        )
    nc.write(
        xyz * 10,
        time,
        np.asarray([al, bl, cl]).T * 10,
        np.asarray([alpha, beta, gamma]).T
    )


def _nc_a_traj(info, gen, gen_fn, chunk):
    out_fn = "{outdir}/{gen}.{outext}".format(gen=gen, **info['cnv2'])
    with XTCTrajectoryFile(gen_fn, 'r') as xtc:
        with NetCDFTrajectoryFile(out_fn, 'w') as nc:
            tot_frames = len(xtc)
            for _ in range(tot_frames // chunk):
                _nc_a_chunk(xtc, nc, chunk)

            if tot_frames % chunk != 0:
                _nc_a_chunk(xtc, nc, chunk=None)

    return out_fn


def _cnv2(info, *, chunk=100):
    info['cnv2'] = {
        'date': datetime.now().isoformat(),
        'chunk': chunk,
        'log': "{workdir}/cnv2.log".format(**info['path']),
        'outdir': "{workdir}/cnv2".format(**info['path']),
        'outext': 'nc',
        'gens': [] if 'cnv1' not in info else info['cnv1']['gens'],
    }

    if not info['cnv1']['success']:
        info['cnv2']['success'] = False
        return info

    log.debug("CNV2: {meta[project]}-{meta[run]}-{meta[clone]}. "
              "Converting to nc"
              .format(**info))

    if info['cnv1']['skip']:
        prev_gens = info['raw']['gens']
    else:
        prev_gens = info['cnv1']['gens']

    done = set(gen for gen, _ in prev_gens)
    assert len(done) == len(prev_gens)
    os.makedirs(info['cnv2']['outdir'], exist_ok=True)
    for gen, gen_fn in prev_gens:
        if gen in done:
            continue
        out_fn = _nc_a_traj(info, gen, gen_fn, chunk)
        info['cnv2']['gens'] += [(gen, out_fn)]

    info['cnv']['success'] = True
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
        topology = "{raw[indir]}/topol.tpr".format(**info),
    else:
        skip = True

    return _cnv1(
        info,
        stride=stride,
        topology=topology,
        skip=skip
    )


def cnv2(info, projcode):
    return _cnv2(info)
