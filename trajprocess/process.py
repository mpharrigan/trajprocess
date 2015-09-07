import subprocess
import os
import glob
import re
from datetime import date

import logging

log = logging.getLogger(__name__)


def nfo_traj(info, *, rncln_re):
    rncln_ma = rncln_re.search(info['raw']['indir'])
    meta = {
        'project': info['meta']['project'],
        'run': int(rncln_ma.group(1)),
        'clone': int(rncln_ma.group(2)),
    }
    path = {'workdir': "processed/{project}/{run}/{clone}".format(**meta)}
    path['info'] = "{workdir}/info.json".format(**path)
    info['meta'] = meta
    info['path'] = path

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


def cat_traj(info, *, gen_glob, gen_re):
    cat_info = {
        'xtc_out': "{workdir}/cat.xtc".format(**info['path']),
        'log_out': "{workdir}/cat.log".format(**info['path']),
    }

    fns = glob.glob(gen_glob.format(**info))
    gen_re = re.compile(gen_re.format(**info))
    gens = sorted(int(gen_re.match(fn).group(1)) for fn in fns)
    cat_info['gen'] = gens[-1] + 1

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

    # Refuse too-short trajectories
    if len(fns) < 2:
        info['cat']['success'] = False
        return info

    # Run trjcat
    with open(info['cat']['log_out'], 'w') as logf:
        subprocess.call(
            ["gmx", "trjcat", "-f"] + fns + ['-o', info['cat']['xtc_out']],
            stdout=logf,
            stderr=subprocess.STDOUT
        )
    info['cat']['success'] = True
    return info


def cat_a4(info):
    return cat_traj(
        info,
        gen_glob="{raw[indir]}/frame*.xtc",
        gen_re="{raw[indir]/frame([0-9]+).xtc",
    )


def cat_21(info):
    return cat_traj(
        info,
        gen_glob="{raw[indir]}/results-???/positions.xtc",
        gen_re="{raw[indir]}/results-([0-9]+)/positions.xtc"
    )


def cnv_traj(info, *, stride=1):
    info['cnv'] = {
        'stride': stride,
        'xtc_out': "{workdir}/cnv.xtc".format(**info['path']),
        'log_out': "{workdir}/cnv.log".format(**info['path']),
    }

    with open(info['cnv']['log_out'], 'w') as logf:
        popen = subprocess.Popen(
            ['gmx', 'trjconv', '-f', info['cat']['xtc_out'], '-o',
             info['cnv']['xtc_out'], '-s',
             '{raw[indir]}/frame0.tpr'.format(**info), '-pbc', 'mol', '-center',
             '-skip', "{cnv[stride]}".format(**info)],
            stdin=subprocess.PIPE,
            stdout=logf,
            stderr=subprocess.STDOUT
        )
        # Center based on 1 - Protein
        # Output 0 - System
        popen.communicate(b"1\n0")

    info['cnv']['success'] = True
    return info


def cnv_21(info):
    info['cnv'] = {
        'stride': 1,
        'xtc_out': "{cat[xtc_out]}".format(**info),
        'success': False,
    }
    return info


def cnv_a4(info):
    if info['meta']['project'] == 'p9752':
        stride = 4
    else:
        stride = 1

    return cnv_traj(
        info,
        stride=stride,
    )
