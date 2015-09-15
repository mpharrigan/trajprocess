inds = range(5)

import subprocess
import os

templ2 = """
parm nocl.nona.nomypc.nowat.{i}.prmtop
trajin nocl.{i}.rst7
autoimage
center @CA
image
trajout done.{i}.pdb
"""


for i in inds:
    with open("cpptraj.tmp", 'w') as f:
        f.write(templ2.format(i=i))
    subprocess.check_call(['cpptraj', '-i', 'cpptraj.tmp'])

os.remove("cpptraj.tmp")
