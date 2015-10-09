"""Functions for performing the individual processing steps.

 - "nfo": Prepare meta-information for each trajectory
 - "cat": Concatenate trajectory parts
 - "cnv": Run trjconv for rough, pbc imaging (gromacs runs only)

"""

import subprocess
import os
import glob
import re
import json
import logging

from mdtraj.formats import XTCTrajectoryFile, NetCDFTrajectoryFile
import mdtraj.utils
import numpy as np

log = logging.getLogger(__name__)


def nfo_traj(info, *, rncln_re, clone=None):
    rncln_ma = rncln_re.search(info['raw']['indir'])
    meta = {
        'project': info['meta']['project'],
        'run': int(rncln_ma.group(1)),
        'clone': int(rncln_ma.group(2)) if clone is None else clone,
    }
    path = {'workdir': "processed/{project}/{run}/{clone}".format(**meta)}
    path['info'] = "{workdir}/info.json".format(**path)
    info['meta'] = meta
    info['path'] = path

    struct_fn = "structs-{meta[project]}.json".format(**info)
    try:
        with open(struct_fn) as f:
            stru = json.load(f)
            string_key = str(info['meta']['run'])  # ugh
            info['top'] = stru[string_key]
    except Exception as e:
        log.warning("No structure information. {}".format(e))

    os.makedirs(path['workdir'], exist_ok=True)
    log.debug("NFO: {project} run {run} clone {clone}".format(**meta))
    return info


def nfo_a4(info):
    return nfo_traj(
        info,
        rncln_re=re.compile(r"RUN(\d+)/CLONE(\d+)/"),
    )


def nfo_21(info):
    return nfo_traj(
        info,
        rncln_re=re.compile(r"RUN(\d+)/CLONE(\d+)/"),
    )


def nfo_bw(info):
    return nfo_traj(
        info,
        rncln_re=re.compile(r"run-(\d+)/"),
        clone=0
    )


def cat_traj(info, *, gen_glob, gen_re):
    cat_info = {
        'xtc_out': "{workdir}/cat.xtc".format(**info['path']),
        'log_out': "{workdir}/cat.log".format(**info['path']),
    }

    fns = glob.glob(gen_glob.format(**info))

    # Refuse too-short trajectories
    if len(fns) < 2:
        if 'cat' not in info:
            info['cat'] = cat_info
        info['cat']['success'] = False
        return info

    gen_re = re.compile(os.path.normpath(gen_re.format(**info)))
    gens = sorted(int(gen_re.match(fn).group(1)) for fn in fns)
    cat_info['gen'] = gens[-1] + 1

    if 'cat' in info and info['cat']['success']:
        # Set up for appending
        conv_fns = [fn for fn in fns
                    if int(gen_re.match(fn).group(1)) >= info['cat']['gen']]
        if len(conv_fns) == 0:
            return info
        conv_fns.append(info['cat']['xtc_out'])
    else:
        conv_fns = fns

    info['cat'] = cat_info

    # Make sure gen indexing matches number of files
    if len(fns) != gens[-1] + 1:
        info['cat']['n_files'] = len(fns)
        log.error("CAT: {meta[project]}-{meta[run]}-{meta[clone]} "
                  "Non contiguous trajectories? "
                  "By regex: {cat[gen]}. Files: {cat[n_files]}".format(**info))
        info['cat']['success'] = False
        return info

    # Give some info
    log.debug("CAT: {meta[project]}-{meta[run]}-{meta[clone]} "
              "found {cat[gen]} trajectories".format(**info))


    # Run trjcat
    with open(info['cat']['log_out'], 'w') as logf:
        subprocess.check_call(
            (["gmx", "trjcat", "-f"]
             + conv_fns + ['-o', info['cat']['xtc_out']]),
            stdout=logf,
            stderr=subprocess.STDOUT
        )
    info['cat']['success'] = True
    return info


def cat_a4(info):
    return cat_traj(
        info,
        gen_glob="{raw[indir]}/frame*.xtc",
        gen_re="{raw[indir]}/frame([0-9]+).xtc",
    )


def cat_21(info):
    return cat_traj(
        info,
        gen_glob="{raw[indir]}/results-???/positions.xtc",
        gen_re="{raw[indir]}/results-([0-9]+)/positions.xtc"
    )


def cat_bw(info):
    info['cat'] = {
        'success': True,
        'gen': 1,
        'xtc_out': "{raw[indir]}/traj_comp.xtc".format(**info)
    }
    return info


def cnv_traj(info, *, stride, topology, skip=False):
    if not info['cat']['success']:
        info['cnv'] = {'success': False}
        log.debug("CNV: {meta[project]}-{meta[run]}-{meta[clone]}. No input"
                  .format(**info))
        return info

    if skip:
        info['cnv'] = {
            'stride': stride,
            'xtc_out': "{cat[xtc_out]}".format(**info),
            'success': True,
        }
        log.debug("CNV: {meta[project]}-{meta[run]}-{meta[clone]}. Skipping"
                  .format(**info))
        return info

    info['cnv'] = {
        'stride': stride,
        'xtc_out': "{workdir}/cnv.xtc".format(**info['path']),
        'log_out': "{workdir}/cnv.log".format(**info['path']),
    }
    log.debug("CNV: {meta[project]}-{meta[run]}-{meta[clone]}. Doing"
              .format(**info))

    with open(info['cnv']['log_out'], 'w') as logf:
        popen = subprocess.Popen(
            ['gmx', 'trjconv', '-f', info['cat']['xtc_out'], '-o',
             info['cnv']['xtc_out'], '-s',
             topology, '-pbc', 'mol', '-center',
             '-skip', "{cnv[stride]}".format(**info)],
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

    info = cnv_to_nc(info)
    info['cnv']['success'] = True
    return info


def _do_a_chunk(xtc, nc, chunk):
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


def cnv_to_nc(info, *, chunk=100):
    log.debug("CNV: {meta[project]}-{meta[run]}-{meta[clone]}. Converting to nc"
              .format(**info))
    info['cnv']['nc_out'] = '{workdir}/cnv.nc'.format(**info['path'])
    with XTCTrajectoryFile(info['cnv']['xtc_out'], 'r') as xtc:
        with NetCDFTrajectoryFile(info['cnv']['nc_out'], 'w') as nc:
            for chunki in range(len(xtc) // chunk):
                _do_a_chunk(xtc, nc, chunk)
            _do_a_chunk(xtc, nc, chunk=None)

    return info


def cnv_21(info):
    return cnv_traj(
        info,
        stride=1,
        topology="",
        skip=True
    )


def cnv_a4(info):
    if info['meta']['project'] == 'p9752':
        stride = 4
    elif info['meta']['project'] == 'p9761':
        stride = 8
    else:
        stride = 1

    return cnv_traj(
        info,
        stride=stride,
        topology="{raw[indir]}/frame0.tpr".format(**info),
    )


def cnv_bw(info):
    return cnv_traj(
        info,
        stride=1,
        topology="{raw[indir]}/topol.tpr".format(**info),
    )
