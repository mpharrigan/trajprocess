import os
from tempfile import mkdtemp
import shutil

from pkg_resources import resource_filename


def mock_project():
    global wd
    wd = mkdtemp()
    os.chdir(wd)

    refdir = resource_filename(__name__, "mock4-reference")
    shutil.copytree(refdir, "trajprocess-testing")
    os.chdir("trajprocess-testing")

    os.mkdir("trajprocess.in.4")
    os.symlink("../PROJ9704", "trajprocess.in.4/p9704")
    os.symlink("../PROJ9752", "trajprocess.in.4/p9752")
    os.symlink("../structs-p9704.json", "trajprocess.in.4/p9704-structs.json")
    os.symlink("../structs-p9752.json", "trajprocess.in.4/p9752-structs.json")
    os.symlink("../tops-p9704", "trajprocess.in.4/p9704-tops")


def cleanup():
    shutil.rmtree(wd)
