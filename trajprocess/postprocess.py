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
from tempfile import TemporaryDirectory

log = logging.getLogger(__name__)


def _norm_cpptraj(cpptraj_selection):
    # Note: replace minus must be first.
    return cpptraj_selection.lower() \
        .replace("-", "-min") \
        .replace("+", "-pl") \
        .replace(":", "mol-") \
        .replace("@", "atm-")


def call_cpptraj_stp(infn, outfn, logfn, *, removes,
                     num_to_keeps, prmtopdir, outtopdir, struct):
    prevs = [None] + [_norm_cpptraj(remove) for remove in removes]

    # Ugh. cpptraj appends names instead of letting you specify the actual
    # out filename. Keep track of these appended names.
    cumprevs = ['.'.join(prevs[1:j][::-1]) for j in range(1, len(prevs) + 1)]

    template = Template("\n".join([
        "{% if prev is none %}",
        "parm {{topdir}}/{{struct}}.prmtop",
        "trajin {{infn}}",
        "{% else %}",
        "parm {{workdir}}/{{cumprev}}.{{struct}}.prmtop",
        "trajin {{workdir}}/{{prev}}.nc",
        "{% endif %}",
        "solvent {{remove}}",
        "closest {{num}} @CA closestout {{workdir}}/{{curr}}.dat outprefix {{workdir}}/{{curr}}",
        "trajout {{workdir}}/{{curr}}.nc",
        ""
    ]))

    varszip = zip(removes, prevs, prevs[1:], cumprevs, num_to_keeps)
    with TemporaryDirectory(dir='/dev/shm', prefix='trajproc-') as td:
        for vars in varszip:
            remove, prev, curr, cumprev, num = vars
            workfile = "{td}/cpptraj.{curr}.tmp"
            workfile = workfile.format(curr=curr, td=td)

            with open(workfile, 'w') as f:
                f.write(template.render(
                        topdir=prmtopdir, struct=struct, infn=infn,
                        workdir=td, cumprev=cumprev, prev=prev, curr=curr,
                        remove=remove, num=num,
                ))

            with open(logfn, 'a') as logf:
                subprocess.check_call(
                        ['cpptraj', '-i', workfile],
                        stderr=subprocess.STDOUT, stdout=logf
                )

        # Move results
        tmp_fn = "{td}/{final}.nc".format(td=td, final=prevs[-1])
        shutil.move(tmp_fn, outfn)

        # Move prmtop if it doesn't already exist
        outtop = ("{outtopdir}/{struct}.strip.prmtop"
                  .format(outtopdir=outtopdir, struct=struct))
        if not os.path.exists(outtop):
            os.makedirs(outtopdir, exist_ok=True)
            tmp_fn = ("{td}/{cumfinal}.{struct}.prmtop"
                      .format(td=td, cumfinal=cumprevs[-1], struct=struct))
            shutil.move(tmp_fn, outtop)


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

    print(removes, num_to_keeps, topdir)


def _call_cpptraj_ctr(info, template, gen, gen_fn):
    workfile = "{outdir}/cpptraj.tmp".format(**info['ctr'])
    out_fn = "{outdir}/{g}.nc".format(g=gen, **info['ctr'])

    log.debug("Centering {}".format(out_fn))

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
    log.info("CTR: {meta[project]}-{meta[run]}-{meta[clone]}. "
             "Using cpptraj to image and center. "
             "Done {done}, todo {todo}"
             .format(done=done, todo=len(info['stp']['gens']) - done, **info))
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
