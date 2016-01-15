from trajprocess.multiprocessing import Pool
from trajprocess.app import execute_task
from trajprocess.systems.nav_no_cpptraj import NaV

with Pool(14) as lbv:
    execute_task(NaV(), lbv)
