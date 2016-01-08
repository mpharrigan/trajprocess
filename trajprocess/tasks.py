import os
import subprocess
import re
import glob

from .prc import PRC
from .project import parse_project, parse_projtype
from .app import config


class Task:
    @property
    def depends(self):
        yield from []

    @property
    def is_done(self):
        return True

    def do(self, tasks):
        return


class RawXTC(Task):
    def __init__(self, prc):
        self.prc = prc

    @property
    def fn(self):
        return "{prc:raw}".format(prc=self.prc)


class PRCTask(Task):
    code = "unsp"
    fext = 'nc'
    dep_class = RawXTC
    needs_log = False

    def __init__(self, prc):
        self.prc = prc

    @property
    def depends(self):
        yield self.dep_class(self.prc)

    @property
    def fn(self):
        return ("{outdir}/{prc:dir}/{code}/{prc:gen}.{fext}"
                .format(outdir=config.outdir,
                        code=self.code,
                        prc=self.prc,
                        fext=self.fext)
                )

    @property
    def log_fn(self):
        return ("{outdir}/{prc:dir}/{code}/{prc:gen}.{fext}"
                .format(outdir=config.outdir,
                        code=self.code,
                        prc=self.prc,
                        fext='log')
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
        if self.needs_log:
            self.do_file(in_fn, out_fn, self.log_fn)
        else:
            self.do_file(in_fn, out_fn)

    def __str__(self):
        return "<{} {}>".format(self.prc, self.code)

    def do_file(self, infn, outfn, logfn=None):
        raise NotImplementedError


class Project(Task):
    dep_class = PRCTask
    projtype = 'x21'

    def __init__(self, project):
        self.project = project
        self.projcode, self.projnum = parse_project(project)
        parse_projtype(self.projtype)

        self.indir = "{indir}/{project}".format(indir=config.indir,
                                                project=project)
        if not os.path.exists(self.indir):
            raise ValueError("Project in directory doesn't exist. "
                             "Looking for {}".format(self.indir))

        self._depends = None

    def _get_prcs(self):
        for run, clone, prc_dir in self.get_run_clones(self.indir):
            for gen, rawfn in self.get_gens(prc_dir):
                yield self._configure(PRC(self.project, run, clone, gen, rawfn))

    def _configure(self, prc):
        return prc

    def get_run_clones(self, indir):
        for fn in glob.iglob("*/"):
            yield 0, 0, fn

    def get_gens(self, prc_dir):
        for fn in glob.iglob("{prc_dir}/*.xtc".format(prc_dir=prc_dir)):
            yield 0, fn

    @property
    def depends(self):
        if self._depends is None:
            self._depends = list(self.dep_class(prc)
                                 for prc in self._get_prcs())
        yield from self._depends


class FahProject(Project):
    prc_glob = "{indir}/RUN*/CLONE*/"
    prc_re = r"{indir}/RUN(\d+)/CLONE(\d+)/"
    gen_re = re.compile("")
    gen_glob = ""

    def get_run_clones(self, indir):
        for fn in glob.iglob(self.prc_glob.format(indir=indir)):
            ma = re.match(self.prc_re.format(indir=indir), fn)
            yield int(ma.group(1)), int(ma.group(2)), fn

    def get_gens(self, prc_dir):
        for fn in (glob.iglob(self.gen_glob.format(prc_dir=prc_dir))):
            yield int(self.gen_re.search(fn).group(1)), fn


class Projectx21(FahProject):
    gen_re = re.compile(r"results-(\d\d\d)/")
    gen_glob = "{prc_dir}/results-???/positions.xtc"


class ProjectxA4(FahProject):
    gen_re = re.compile(r"frame(\d+).xtc")
    gen_glob = "{prc_dir}/frame*.xtc"

    def _configure(self, prc):
        prc.meta['needs_trjconv'] = True
        prc.meta['has_overlapping_frames'] = True
        prc.meta['tpr_fn'] = ("{indir}/frame0.tpr"
                              .format(indir=os.path.dirname(prc.in_fn)))
        return prc


class ProjectBluewaters(Project):
    def _configure(self, prc):
        prc.meta['needs_trjconv'] = True
        prc.meta['tpr_fn'] = ("{indir}/topol.tpr"
                              .format(indir=os.path.dirname(prc.in_fn)))
        return prc

    def get_run_clones(self, indir):
        # TODO
        for fn in glob.iglob("*/"):
            yield 0, 0, fn

    def get_gens(self, prc_dir):
        # TODO
        for fn in glob.iglob("{prc_dir}/*.xtc".format(prc_dir=prc_dir)):
            yield 0, fn
