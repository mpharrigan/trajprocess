#!/usr/bin/env python
import pprint as pp
from ipyparallel import Client

c = Client()
qs = c.queue_status()
pp.pprint(qs)

completed = 0
in_prog = 0

for k in qs:
    try:
        k = int(k)
        val = qs[k]
        completed += val['completed']
        in_prog += val['queue'] + val['tasks']
    except ValueError:
        pass

total = completed + in_prog + qs['unassigned']

for lab, val in zip(
        ["Completed", "In Progress", "Total"],
        [completed, in_prog, total]):
    print("{:20s} {:7}".format(lab, val))
