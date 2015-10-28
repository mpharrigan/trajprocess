"""Perform additional, optional processing steps on trajectories

 - "stp": strip all but closest solvent molecules
 - "ctr": run autoimage and center

"""

import subprocess
import os
import shutil
import logging
from datetime import datetime

from jinja2 import Template

log = logging.getLogger(__name__)


def _norm_cpptraj(cpptraj_selection):
    # Note: replace minus must be first.
    return cpptraj_selection.lower() \
        .replace("-", "-min") \
        .replace("+", "-pl") \
        .replace(":", "mol-") \
        .replace("@", "atm-")


def _call_cpptraj_stp(info, template, gen, removes, prevs, cumprevs,
                      num_to_keeps):
    # Warning: I use workdir as the name for the proj-run-clone directory too
    # sorry
    workdir = "{outdir}/{g}".format(g=gen, **info['stp'])
    os.makedirs(workdir, exist_ok=True)

    varszip = zip(removes, prevs, prevs[1:], cumprevs, num_to_keeps)
    for vars in varszip:
        remove, prev, curr, cumprev, num = vars
        workfile = "{workdir}/cpptraj.{curr}.tmp"
        workfile = workfile.format(curr=curr, workdir=workdir)

        with open(workfile, 'w') as f:
            f.write(template.render(
                remove=remove, prev=prev, curr=curr, cumprev=cumprev, num=num,
                workdir=workdir, g=gen, **info
            ))

        with open(info['stp']['log'], 'a') as logf:
            subprocess.check_call(
                ['cpptraj', '-i', workfile],
                stderr=subprocess.STDOUT, stdout=logf
            )

    # Move results
    tmp_fn = "{workdir}/{final}.nc".format(workdir=workdir, final=prevs[-1])
    out_fn = "{outdir}/{g}.nc".format(g=gen, **info['stp'])
    shutil.move(tmp_fn, out_fn)

    # Move prmtop if it doesn't already exist
    if not os.path.exists(info['stp']['outtop']):
        tmp_fn = ("{workdir}/{cumfinal}.{struct}.prmtop"
                  .format(workdir=workdir, cumfinal=cumprevs[-1],
                          struct=info['top']['struct']))
        shutil.move(tmp_fn, info['stp']['outtop'])


    # Remove files
    shutil.rmtree(workdir)
    return out_fn


def _stp(info, *, removes, num_to_keeps, topdir):
    info['stp'] = {
        'log': "{workdir}/stp.log".format(**info['path']),
        'outdir': "{workdir}/stp".format(**info['path']),
        'topdir': topdir,
        'removes': removes,
        'num_to_keeps': num_to_keeps,
        'date': datetime.now().isoformat(),
        'gens': [] if 'stp' not in info else info['stp']['gens'],
    }
    info['stp']['outtop'] = ("{stp[topdir]}/{top[struct]}.strip.prmtop"
                             .format(**info))

    if not info['cnv2']['success']:
        info['stp']['success'] = False
        return info

    log.debug("STP: {meta[project]}-{meta[run]}-{meta[clone]}. Doing"
              .format(**info))
    prevs = [None] + [_norm_cpptraj(remove) for remove in removes]

    # Ugh. cpptraj appends names instead of letting you specify the actual
    # out filename. Keep track of these appended names.
    cumprevs = ['.'.join(prevs[1:j][::-1]) for j in range(1, len(prevs) + 1)]

    template = Template("\n".join([
        "{% if prev is none %}",
        "parm {{stp['topdir']}}/{{top['struct']}}.prmtop",
        "trajin {{cnv2['gens'][g]}}",
        "{% else %}",
        "parm {{workdir}}/{{cumprev}}.{{top['struct']}}.prmtop",
        "trajin {{workdir}}/{{prev}}.nc",
        "{% endif %}",
        "solvent {{remove}}",
        "closest {{num}} @CA closestout {{workdir}}/{{curr}}.dat outprefix {{workdir}}/{{curr}}",
        "trajout {{workdir}}/{{curr}}.nc",
        ""
    ]))

    done = len(info['stp']['gens'])
    os.makedirs(info['stp']['outdir'], exist_ok=True)
    for gen, gen_fn in enumerate(info['cnv2']['gens']):
        if gen < done:
            continue
        out_fn = _call_cpptraj_stp(info, template, gen, removes, prevs,
                                   cumprevs, num_to_keeps)
        info['stp']['gens'] += [out_fn]

    info['stp']['success'] = True

    return info


def stp(info, systemcode):
    if systemcode == 'nav':

        removes = [":WAT", ":MY", "@Na+", "@Cl-"]
        num_to_keeps = [10000, 100, 20, 20]
        topdir = "tops-p9704"
    elif systemcode == 'trek':
        removes = [":WAT", ":PC", ":PE", "@K+", "@Cl-"]
        num_to_keeps = [5000, 30, 30, 20, 20]
        topdir = "tops-p9712"
    else:
        raise ValueError

    return _stp(
        info,
        removes=removes,
        num_to_keeps=num_to_keeps,
        topdir=topdir,
    )


def _call_cpptraj_ctr(info, template, gen, gen_fn):
    workfile = "{outdir}/cpptraj.tmp".format(**info['ctr'])
    out_fn = "{outdir}/{g}.nc".format(g=gen, **info['ctr'])

    with open(workfile, 'w') as f:
        f.write(template.format(gen_fn=gen_fn, out_fn=out_fn, **info))

    with open(info['ctr']['log'], 'a') as logf:
        subprocess.check_call(
            ['cpptraj', '-i', workfile],
            stderr=subprocess.STDOUT, stdout=logf
        )
    os.remove(workfile)
    return out_fn


def _ctr(info):
    info['ctr'] = {
        'log': "{workdir}/ctr.log".format(**info['path']),
        'outdir': "{workdir}/ctr".format(**info['path']),
        'date': datetime.now().isoformat(),
        'gens': [] if 'ctr' not in info else info['ctr']['gens'],
    }

    if not info['stp']['success']:
        info['ctr']['success'] = False
        return info

    log.debug("CTR: {meta[project]}-{meta[run]}-{meta[clone]}. Doing"
              .format(**info))

    template = "\n".join([
        "parm {stp[outtop]}",
        "trajin {gen_fn}",
        "autoimage",
        "center @CA",
        "image",
        "trajout {out_fn}",
        "",
    ])

    done = len(info['ctr']['gens'])
    os.makedirs(info['ctr']['outdir'], exist_ok=True)
    for gen, gen_fn in enumerate(info['stp']['gens']):
        if gen < done:
            continue
        out_fn = _call_cpptraj_ctr(info, template, gen, gen_fn)
        info['ctr']['gens'] += [out_fn]

    info['ctr']['success'] = True
    return info


def ctr(info, systemcode):
    return _ctr(info)
