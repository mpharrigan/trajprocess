import os
import subprocess

from .project import parse_project, parse_projtype, get_gens, get_prcs
from .app import config


class Task:
    @property
    def is_done(self):
        return False

    def do(self, tasks):
        return


class PRCTask(Task):
    code = "unsp"
    fext = 'nc'

    def __init__(self, prc):
        self.prc = prc

    @property
    def fn(self):
        return ("{outdir}/{prc:dir}/{code}/{prc:gen}.{fext}"
                .format(outdir=config.outdir,
                        code=self.code,
                        prc=self.prc,
                        fext=self.fext)
                )

    @property
    def is_done(self):
        return os.path.exists(self.fn)

    def do(self, task):
        in_fn = task.fn
        if not os.path.exists(in_fn):
            raise FileNotFoundError("Input file for task {} not found. "
                                    "Looking for {}"
                                    .format(self, in_fn))

        out_fn = self.fn
        out_dir = os.path.dirname(out_fn)
        os.makedirs(out_dir, exist_ok=True)
        self.do_file(in_fn, out_fn)

    def __str__(self):
        return "<{} {}>".format(self.prc, self.code)

    def do_file(self, infn, outfn):
        raise NotImplementedError


class RawXTC(Task):
    def __init__(self, prc):
        self.prc = prc

    @property
    def is_done(self):
        return True

    @property
    def fn(self):
        return "{prc:raw}".format(prc=self.prc)

    def do(self, tasks):
        return


class Trjconv(PRCTask):
    code = 'cnv1'
    fext = 'xtc'

    @property
    def depends(self):
        yield RawXTC(self.prc)

    def do_file(self, infn, outfn):
        assert os.path.exists(infn)
        subprocess.call(['touch', outfn])


class ConvertToNC(PRCTask):
    code = 'cnv2'

    @property
    def depends(self):
        if "needs_trjconv" in self.prc.flags:
            yield Trjconv(self.prc)
        else:
            yield RawXTC(self.prc)

    def do_file(self, infn, outfn):
        assert os.path.exists(infn)
        subprocess.call(['touch', outfn])


class Strip(PRCTask):
    code = 'stp'

    @property
    def depends(self):
        yield ConvertToNC(self.prc)

    def do_file(self, infn, outfn):
        assert os.path.exists(infn)
        subprocess.call(['touch', outfn])


class Center(PRCTask):
    code = 'ctr'

    @property
    def depends(self):
        yield Strip(self.prc)

    def do_file(self, infn, outfn):
        assert os.path.exists(infn)
        subprocess.call(['touch', outfn])


class Project(Task):
    def __init__(self, project, projtype):
        self.project = project
        self.projcode, self.projnum = parse_project(project)
        self.projtype = parse_projtype(projtype)

        self.indir = "{indir}/{project}".format(indir=config.indir,
                                                project=project)
        if not os.path.exists(self.indir):
            raise ValueError("Project in directory doesn't exist. "
                             "Looking for {}".format(self.indir))

        self._depends = None

    def _get_depends(self):
        for prc in get_prcs(self.project, self.projtype, self.indir):
            yield Center(prc)

    @property
    def depends(self):
        if self._depends is None:
            self._depends = list(self._get_depends())

        yield from self._depends


class NaV(Task):
    @property
    def depends(self):
        yield from [
            Project("p9704", 'x21'),
            Project("p9752", 'xa4'),
        ]
