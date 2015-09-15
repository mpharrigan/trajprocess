inds = range(5)
remove = [":WAT", ":MY", "@Na+", "@Cl-"]
prevs = ['', 'nowat', 'nomypc',  'nona', 'nocl']
cumprev = ['.'.join(prevs[1:j][::-1]) for j in range(1,len(prevs))]
nums = [10000, 100, 20, 20]

import subprocess
import os
from jinja2 import Template
import sys
import time

#print(cumprev)
#sys.exit(1)


templ = Template("""
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



for i in inds:
    for rem, prev, cprev, cur, num in zip(remove, prevs, cumprev, prevs[1:], nums):
        with open("cpptraj.tmp", 'w') as f:
            f.write(templ.render(i=i, rem=rem, prev=prev, cprev=cprev, cur=cur, num=num))
        subprocess.check_call(['cpptraj', '-i', 'cpptraj.tmp'])
        #print(i, rem, prev, cprev, cur, num)
        #time.sleep(1)

os.remove("cpptraj.tmp")
