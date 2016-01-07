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

    def __init__(self, prc, gen):
        self.prc = prc
        self.gen = gen

    @property
    def fn(self):
        return ("{outdir}/{prc:dir}/{code}/{gen}.{fext}"
                .format(outdir=config.outdir,
                        code=self.code,
                        prc=self.prc,
                        gen=self.gen,
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
        return "<{} {} {}>".format(self.prc, self.gen, self.code)

    def do_file(self, infn, outfn):
        raise NotImplementedError


class RawXTC(Task):
    @property
    def is_done(self):
        return True

    @property
    def fn(self):
        return "{indir}/RUN{run}/CLONE{clone}/frame{gen}.xtc"

    def do(self, tasks):
        return


class Trjconv(PRCTask):
    code = 'cnv1'
    fext = 'xtc'

    @property
    def depends(self):
        yield RawXTC()

    def do_file(self, infn, outfn):
        assert os.path.exists(infn)
        subprocess.call(['touch', outfn])


class ConvertToNC(PRCTask):
    code = 'cnv2'

    @property
    def depends(self):
        self.is_gromacs = True
        if self.is_gromacs:
            yield Trjconv(self.prc, self.gen)
        else:
            yield RawXTC()

    def do_file(self, infn, outfn):
        assert os.path.exists(infn)
        subprocess.call(['touch', outfn])


class Strip(PRCTask):
    code = 'stp'

    @property
    def depends(self):
        yield ConvertToNC(self.prc, self.gen)

    def do_file(self, infn, outfn):
        assert os.path.exists(infn)
        subprocess.call(['touch', outfn])


class Center(PRCTask):
    code = 'ctr'

    @property
    def depends(self):
        yield Strip(self.prc, self.gen)

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
            for gen in get_gens(prc, self.projtype):
                yield Center(prc, gen)

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
