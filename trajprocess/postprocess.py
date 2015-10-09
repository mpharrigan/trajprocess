"""Perform additional, optional processing steps on trajectories

 - "stp": strip all but closest solvent molecules
 - "ctr": run autoimage and center

"""

import subprocess
import os
import shutil
import glob

from jinja2 import Template


def _norm_cpptraj(cpptraj_selection):
    # Note: replace minus must be first.
    return cpptraj_selection.lower() \
        .replace("-", "-min") \
        .replace("+", "-pl") \
        .replace(":", "mol-") \
        .replace("@", "atm-")


def stp_traj(info, *, removes, num_to_keeps, topdir, topext="prmtop",
             trajext='nc'):
    if not info['cnv']['success']:
        info['stp'] = {'success': False}
        return info

    prevs = [None] + [_norm_cpptraj(remove) for remove in removes]

    # Ugh. cpptraj appends names instead of letting you specify the actual
    # out filename. Keep track of these appended names.
    cumprevs = ['.'.join(prevs[1:j][::-1]) for j in range(1, len(prevs) + 1)]

    template = Template("\n".join([
        "{% if prev is none %}",
        "parm {{topdir}}/{{top['struct']}}.{{topext}}",
        "trajin {{cnv['nc_out']}}",
        "{% else %}",
        "parm {{stp['cpp_workdir']}}/{{cumprev}}.{{top['struct']}}.prmtop",
        "trajin {{stp['cpp_workdir']}}/{{prev}}.{{trajext}}",
        "{% endif %}",
        "solvent {{remove}}",
        "closest {{num}} @CA closestout {{stp['cpp_workdir']}}/{{curr}}.dat outprefix {{stp['cpp_workdir']}}/{{curr}}",
        "trajout {{stp['cpp_workdir']}}/{{curr}}.{{trajext}}",
        ""
    ]))

    info['stp'] = {
        'cpp_workdir': "{path[workdir]}/cpptraj".format(**info),
        'cpp_archive': "{path[workdir]}/cpptraj.tar.gz".format(**info),
        'log_out': "{path[workdir]}/stp.log".format(**info),
        'removes': removes,
        'num_to_keeps': num_to_keeps,
    }

    os.makedirs(info['stp']['cpp_workdir'], exist_ok=True)

    varszip = zip(removes, prevs, prevs[1:], cumprevs, num_to_keeps)
    for vars in varszip:
        remove, prev, curr, cumprev, num = vars
        workfile = "{stp[cpp_workdir]}/cpptraj.{curr}.tmp"
        workfile = workfile.format(curr=curr, **info)

        with open(workfile, 'w') as f:
            f.write(template.render(
                remove=remove, prev=prev, curr=curr, cumprev=cumprev, num=num,
                topdir=topdir, topext=topext, trajext=trajext, **info
            ))

        with open(info['stp']['log_out'], 'a') as logf:
            subprocess.check_call(
                ['cpptraj', '-i', workfile],
                stderr=subprocess.STDOUT, stdout=logf
            )

    # Move results
    assert trajext == 'nc'
    info['stp']['nc_out'] = "{path[workdir]}/stp.nc".format(**info)
    info['stp']['prmtop'] = "{path[workdir]}/stp.prmtop".format(**info)
    shutil.move(
        "{stp[cpp_workdir]}/{final}.{trajext}"
            .format(final=prevs[-1], trajext=trajext, **info),
        "{stp[nc_out]}".format(**info)
    )
    shutil.move(
        "{stp[cpp_workdir]}/{cumfinal}.{top[struct]}.prmtop"
            .format(cumfinal=cumprevs[-1], **info),
        "{stp[prmtop]}".format(**info)
    )

    # Remove intermediate trajectory files
    for fn in glob.glob("{stp[cpp_workdir]}/*.nc".format(**info)):
        os.remove(fn)

    # Tar up workdir
    subprocess.check_call(
        ['tar', '-czf', info['stp']['cpp_archive'], info['stp']['cpp_workdir']],
        stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL
    )

    # Clean up workdir
    shutil.rmtree(info['stp']['cpp_workdir'])

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
        removes=[":WAT", ":PC", ":PE", "@K+", "@Cl-"],
        num_to_keeps=[5000, 30, 30, 20, 20],
        topdir="tops-p9712",
    )


def ctr_traj(info):
    if not info['stp']['success']:
        info['ctr'] = {'success': False}
        return info

    info['ctr'] = {
        'nc_out': "{workdir}/ctr.nc".format(**info['path']),
        'log_out': "{workdir}/ctr.log".format(**info['path']),
        'prmtop': info['stp']['prmtop'],
    }

    template = "\n".join([
        "parm {stp[prmtop]}",
        "trajin {stp[nc_out]}",
        "autoimage",
        "center @CA",
        "image",
        "trajout {ctr[nc_out]}",
        "",
    ])

    workfile = "{workdir}/cpptraj.tmp".format(**info['path'])
    with open(workfile, 'w') as f:
        f.write(template.format(**info))

    with open(info['ctr']['log_out'], 'w') as logf:
        subprocess.check_call(
            ['cpptraj', '-i', workfile],
            stderr=subprocess.STDOUT, stdout=logf
        )

    os.remove(workfile)
    return info
