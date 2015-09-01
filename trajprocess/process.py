import subprocess
import os
import glob
import re
from datetime import date


def nfo_traj(info, *, rncln_re):
    rncln_ma = rncln_re.search(info['raw_indir'])
    info.update({
        'nfo_indir': os.path.realpath(info['raw_indir']),
        'run': int(rncln_ma.group(1)),
        'clone': int(rncln_ma.group(2)),
        'idate': date.today().isoformat(),
    })
    info['nfo_outdir'] = "{idate}/{project}/{run}/{clone}/".format(**info)
    info['nfo_nfoout'] = "{nfo_outdir}".format(**info)
    os.makedirs(info['nfo_outdir'], exist_ok=True)
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


def cat_traj(info, *, gen_glob):
    info['cat_outdir'] = "{nfo_outdir}".format(**info)
    info['cat_xtcout'] = "{cat_outdir}/cat.xtc".format(**info)
    info['cat_logout'] = "{cat_outdir}/cat.log".format(**info)
    os.makedirs(info['cat_outdir'], exist_ok=True)

    fns = glob.glob(gen_glob.format(**info))

    # Refuse too-short trajectories
    if len(fns) < 2:
        info['cat_success'] = False
        return info

    # Run trjcat
    with open(info['cat_logout'], 'w') as logf:
        subprocess.call(
            ["gmx", "trjcat", "-f"] + fns + ['-o', info['cat_xtcout']],
            stdout=logf,
            stderr=subprocess.STDOUT
        )
    info['cat_success'] = True
    return info


def cat_a4(info):
    return cat_traj(
        info,
        gen_glob="{cat_indir}/frame*.xtc",
    )


def cat_21(info):
    return cat_traj(
        info,
        gen_glob="{cat_indir}/results-???/positions.xtc",
    )


def cnv_traj(info, *, stride=1, do=True):
    info['stride'] = stride
    info['cnv_outdir'] = "{cat_outdir}".format(**info)
    info['cnv_xtcout'] = "{cnv_outdir}/cnv.xtc".format(**info)
    info['cnv_logout'] = "{cnv_outdir}/cnv.log".format(**info)
    os.makedirs(info['cnv_outdir'], exist_ok=True)

    if not do:
        info['cnv_success'] = False
        return info

    with open(info['cnv_logout'], 'w') as logf:
        popen = subprocess.Popen(
            ['gmx', 'trjconv', '-f', info['cat_xtcout'], '-o',
             info['cnv_xtcout'], '-s', '{cat_indir}/frame0.tpr'.format(**info),
             '-pbc', 'mol', '-center', '-skip', "{stride}".format(**info)],
            stdin=subprocess.PIPE,
            stdout=logf,
            stderr=subprocess.STDOUT
        )
        # Center based on 1 - Protein
        # Output 0 - System
        popen.communicate(b"1\n0")

    info['cnv_success'] = True
    return info


def cnv_21(info):
    return cnv_traj(
        info,
        stride=1,
        do=False,
    )


def cnv_a4(info):
    if info['project'] == 'p9752':
        stride = 4
    else:
        stride = 1

    return cnv_traj(
        info,
        stride=stride,
    )
