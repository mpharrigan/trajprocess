from pathlib import Path

from nose import with_setup

from .mock3 import mock_project, cleanup
import trajprocess.files
import trajprocess.cleanup


@with_setup(mock_project, cleanup)
def test_cleanup():
    # Run project and clean up
    trajprocess.files.main_trek()
    trajprocess.cleanup.cleanup()

    home = Path("processed.v2")

    # Make sure we deleted stuff
    print("\n".join(str(s) for s in Path("processed.v2").glob("**/*")))
    assert len(list(home.glob("*/*/*/cnv1/*"))) == 0
    assert len(list(home.glob("*/*/*/cnv2/*"))) == 0
    assert len(list(home.glob("*/*/*/stp/*"))) == 0

    # Make sure cleaning up a clean directory doesn't do bad things
    trajprocess.cleanup.cleanup()

    # Make sure we didn't delete good stuff
    for i in range(3):
        assert (home / "p9712/5/32/ctr/{}.nc".format(i)).exists()
