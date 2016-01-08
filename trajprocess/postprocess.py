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
from tempfile import TemporaryDirectory, NamedTemporaryFile

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


def call_cpptraj_ctr(infn, outfn, logfn, *, stptopdir, struct):
    # Make sure this is sync-ed with above
    topfn = ("{stptopdir}/{struct}.strip.prmtop"
             .format(stptopdir=stptopdir, struct=struct))

    template = "\n".join([
        "parm {topfn}",
        "trajin {infn}",
        "autoimage",
        "center @CA",
        "image",
        "trajout {outfn}",
        "",
    ])

    with NamedTemporaryFile('w', delete=False) as tf:
        tf.write(template.format(topfn=topfn, infn=infn, outfn=outfn))
        tf.flush()
        os.fsync(tf.fileno())

        with open(logfn, 'a') as logf:
            subprocess.check_call(
                    ['cpptraj', '-i', tf.name],
                    stderr=subprocess.STDOUT, stdout=logf
            )
