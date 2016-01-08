from nose import with_setup
from .mock4 import mock_project as mock4, cleanup as cleanup4

from trajprocess.app import execute_task

PROJ9704_FRAMES_PER_GEN = 8
PROJ9752_FRAMES_PER_GEN = 1
N_TRIM_ATOMS = 56700


# TODO: cleanup
@with_setup(mock4, None)
def test_task():
    pass
