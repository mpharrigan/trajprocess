from trajprocess.postprocess import _norm_cpptraj
import os

from nose import with_setup
import mdtraj
from .utils import generate_project, cleanup, generate_bw
from trajprocess.files import Project, record


def test_norm_cpptraj():
    assert _norm_cpptraj(":WAT") == 'mol-wat'
    assert _norm_cpptraj(":MY") == 'mol-my'
    assert _norm_cpptraj("@Na+") == "atm-na-pl"
    assert _norm_cpptraj("@Cl-") == "atm-cl-min"


@with_setup(generate_project, cleanup)
def test_stp():
    project = Project("p1234", "data/PROJ1234", 'xa4')
    raw_infos = list(project.get_infos())
    nfo_infos = list(map(record(project.nfo), raw_infos))
    cat_infos = list(map(record(project.cat), nfo_infos))
    cnv_infos = list(map(record(project.cnv), cat_infos))
    stp_infos = list(map(record(project.stp), cnv_infos))

    assert len(stp_infos) > 0

    for info in stp_infos:
        assert os.path.exists(info['stp']['xtc_out'])

        with mdtraj.open(info['stp']['xtc_out']) as tfile:
            xyz, time, step, box = tfile.read()
            print("Shape", xyz.shape)
            assert xyz.shape == (21, 22, 3), xyz.shape
