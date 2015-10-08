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


def stp_traj(info, *, removes, num_to_removes, topdir, topext="prmtop", trajext='dcd'):
    prevs = [None] + [_norm_cpptraj(remove) for remove in removes]

    # Ugh. cpptraj appends names instead of letting you specify the actual
    # out filename. Keep track of these appended names.
    cumprevs = ['.'.join(prevs[1:j][::-1]) for j in range(1, len(prevs))]

    template = Template("""
    {% if prev is None %}
    parm {{topdir}}/{{top['struct']}}.{{topext}}
    trajin {{cnv['xtc_out']}}
    {% else %}
    parm {{path['workdir']}}/{{cumprev}}.prmtop
    trajin {{path['workdir']}}/{{prev}}.{{trajext}}
    {% endif %}
    solvent {{remove}}
    closest {{num}} @CA closestout {{path['workdir']}}/{{curr}}.dat outprefix {{curr}}
    trajout {{path['workdir']}}/{{curr}}.{{trajext}}
    """)

    varszip = zip(removes, prevs, prevs[1:], cumprevs, num_to_removes)
    for vars in varszip:
        remove, prev, curr, cumprev, num = vars
        tmpfile = "{path[workdir]}/cpptraj.{curr}.tmp".format(curr=curr, **info)
        with open(tmpfile, 'w') as f:
            f.write(template.render(
                remove=remove, prev=prev, curr=curr, cumprev=cumprev, num=num,
                topdir=topdir, topext=topext, trajext=trajext, **info
            ))
        subprocess.check_call(
            ['cpptraj', '-i', 'cpptraj.tmp'],
            cwd=info['path']['workdir']
        )

    # TODO
    return info


def stp_nav(info):
    return stp_traj(
        info,
        remove=[":WAT", ":MY", "@Na+", "@Cl-"],
        num_to_remove=[10000, 100, 20, 20],
    )
