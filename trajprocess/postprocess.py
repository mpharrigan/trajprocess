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


def stp_traj(info, *, remove, num_to_remove):
    prevs = [""] + [_norm_cpptraj(sel) for sel in remove]

    # Ugh. cpptraj appends names instead of letting you specify the actual
    # out filename. Keep track of these appended names.
    cumprev = ['.'.join(prevs[1:j][::-1]) for j in range(1, len(prevs))]

    template = Template("""
    {% if prev == '' %}
    parm {{i}}.prmtop
    trajin {{i}}.rst7
    {% else %}
    parm {{cprev}}.{{i}}.prmtop
    trajin {{prev}}.{{i}}.rst7
    {% endif %}
    solvent {{rem}}
    closest {{num}} @CA closestout {{cur}}.{{i}}.dat outprefix {{cur}}
    trajout {{cur}}.{{i}}.rst7
    """)
    for i in range(len(prevs)):
        for rem, prev, cprev, cur, num in zip(remove, prevs, cumprev,
                                              prevs[1:], num_to_remove):
            with open("cpptraj.tmp", 'w') as f:
                f.write(template.render(i=i, rem=rem,
                                        prev=prev, cprev=cprev,
                                        cur=cur, num=num))
            subprocess.check_call(['cpptraj', '-i', 'cpptraj.tmp'])

    # TODO
    return info


def stp_nav(info):
    return stp_traj(
        info,
        remove=[":WAT", ":MY", "@Na+", "@Cl-"],
        num_to_remove=[10000, 100, 20, 20],
    )
