import os
from tempfile import mkdtemp
from pkg_resources import resource_filename
import shutil


def mock_project():
    global wd
    wd = mkdtemp()
    os.chdir(wd)

    refdir = resource_filename(__name__, "mock2-reference")
    shutil.copytree(refdir, "trajprocess-testing")
    os.chdir("trajprocess-testing")


def cleanup():
    shutil.rmtree(wd)
