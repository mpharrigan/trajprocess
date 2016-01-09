from ipyparallel import Client
from trajprocess.app import execute_task
from trajprocess.systems.nav import NaV

c = Client()
lbv = c.load_balanced_view()
execute_task(NaV(), lbv)
