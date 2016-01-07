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
    def is_done(self):
        return os.path.exists("{outdir}/{prc:dir}/{code}/{gen}.{fext}"
                              .format(outdir=config.outdir,
                                      prc=self.prc,
                                      code=self.code,
                                      gen=self.gen,
                                      fext=self.fext))

    def do(self, prior):
        out_fn = ("{outdir}/{prc:dir}/{code}/{gen}.nc"
                  .format(outdir=config.outdir,
                          code=self.code,
                          prc=self.prc,
                          gen=self.gen)
                  )
        out_dir = os.path.dirname(out_fn)
        os.makedirs(out_dir, exist_ok=True)
        self.do_file(prior, out_fn)

    def do_file(self, infn, outfn):
        raise NotImplementedError


class Strip(Task):
    def __init__(self, prc, gen):
        self.prc = prc
        self.gen = gen

    @property
    def depends(self):
        yield from []

    @property
    def is_done(self):
        return True


class Center(PRCTask):
    code = 'ctr'

    @property
    def depends(self):
        yield from [Strip(self.prc, self.gen)]

    def do_file(self, infn, outfn):
        subprocess.call(['touch', 'outfn'])


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
