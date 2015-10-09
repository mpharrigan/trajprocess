"""Perform additional, optional processing steps on trajectories

 - "stp": strip all but closest solvent molecules and run autoimage

"""

import subprocess
import os
from jinja2 import Template
import sys
import time


def _norm_cpptraj(cpptraj_selection):
    # Note: replace minus must be first.
    return cpptraj_selection.lower() \
        .replace("-", "-min") \
        .replace("+", "-pl") \
        .replace(":", "mol-") \
        .replace("@", "atm-")


def stp_traj(info, *, removes, num_to_keeps, topdir, topext="prmtop",
             trajext='nc'):
    prevs = [None] + [_norm_cpptraj(remove) for remove in removes]

    # Ugh. cpptraj appends names instead of letting you specify the actual
    # out filename. Keep track of these appended names.
    cumprevs = ['.'.join(prevs[1:j][::-1]) for j in range(1, len(prevs))]

    template = Template("\n".join([
        "{% if prev is none %}",
        "parm {{topdir}}/{{top['struct']}}.{{topext}}",
        "trajin {{cnv['nc_out']}}",
        "{% else %}",
        "parm {{path['workdir']}}/{{cumprev}}.prmtop",
        "trajin {{path['workdir']}}/{{prev}}.{{trajext}}",
        "{% endif %}",
        "solvent {{remove}}",
        "closest {{num}} @CA closestout {{path['workdir']}}/{{curr}}.dat outprefix {{curr}}",
        "trajout {{path['workdir']}}/{{curr}}.{{trajext}}",
        ""
    ]))

    varszip = zip(removes, prevs, prevs[1:], cumprevs, num_to_keeps)
    for vars in varszip:
        remove, prev, curr, cumprev, num = vars
        workfile = "{path[workdir]}/cpptraj.{curr}.tmp"
        workfile = workfile.format(curr=curr, **info)

        with open(workfile, 'w') as f:
            f.write(template.render(
                remove=remove, prev=prev, curr=curr, cumprev=cumprev, num=num,
                topdir=topdir, topext=topext, trajext=trajext, **info
            ))
        subprocess.check_call(
            ['cpptraj', '-i', workfile],
            # TODO: stdout, stderr redirect
        )

    # TODO
    return info


def stp_nav(info):
    return stp_traj(
        info,
        removes=[":WAT", ":MY", "@Na+", "@Cl-"],
        num_to_keeps=[10000, 100, 20, 20],
        topdir="tops-p9704",
    )


def stp_trek(info):
    return stp_traj(
        info,
        removes=[":WAT", "@K+"],
        num_to_keeps=[10000, 20],
        topdir="tops-p9712",
    )
