import errno
import glob
import json
import os
import re

from .config import config
from .prc import PRCG
from .project import parse_project


class Task:
    is_ephemeral = False
    is_dummy = False

    @property
    def depends(self):
        yield from []

    @property
    def is_done(self):
        raise NotImplementedError

    def do(self, tasks):
        raise NotImplementedError

    def __str__(self):
        return self.__class__.__name__


class Dummy:
    is_dummy = True

    def do(self, tasks):
        return

    @property
    def is_done(self):
        return True


class RawXTC(Dummy, Task):
    def __init__(self, prcg):
        self.prcg = prcg

    @property
    def fn(self):
        return "{prcg:raw}".format(prcg=self.prcg)

    def __str__(self):
        return "<{} {}>".format(self.__class__.__name__, self.prcg)


class Clean(Task):
    ephemeral_task_classes = []
    dep_class = RawXTC
    delete_logs = False
    delete_empty_dirs = True

    def __init__(self, prcg):
        self.prcg = prcg

    @property
    def depends(self):
        yield self.dep_class(self.prcg)

    @property
    def ephemeral_tasks(self):
        for tsk_class in self.ephemeral_task_classes:
            yield tsk_class(self.prcg)

    @property
    def is_done(self):

        # The following check is important:
        # If none of them have run, then there's nothing to clean up (yet)
        # and it thinks it's done! It will not be scheduled by an async
        # task scheduler.
        immediate_dep = self.dep_class(self.prcg)
        if not immediate_dep.is_done:
            return False

        for task in self.ephemeral_tasks:
            if os.path.exists(task.fn):
                return False

            if self.delete_logs and os.path.exists(task.log_fn):
                return False

        if self.delete_logs and os.path.exists(immediate_dep.log_fn):
            return False

        return True

    def _delete_empty_dirs(self):
        for root, dirs, files in os.walk("{prcg:dir}".format(prcg=self.prcg),
                                         topdown=False):
            for d in dirs:
                try:
                    os.rmdir(os.path.join(root, d))
                except OSError as e:
                    if e.errno == errno.ENOTEMPTY:
                        pass
                    else:
                        print(root, dirs, files, d)
                        raise

    def do(self, tasks):
        for task in self.ephemeral_tasks:
            try:
                os.remove(task.fn)
            except FileNotFoundError:
                pass

            if self.delete_logs:
                try:
                    os.remove(task.log_fn)
                except FileNotFoundError:
                    pass

        if self.delete_logs:
            immediate_dep = self.dep_class(self.prcg)
            try:
                os.remove(immediate_dep.log_fn)
            except FileNotFoundError:
                pass

        if self.delete_empty_dirs:
            self._delete_empty_dirs()

    def __str__(self):
        return "<{} {}>".format(self.__class__.__name__, self.prcg)


class PRCGTask(Task):
    code = "unsp"
    fext = 'nc'
    dep_class = RawXTC
    needs_log = False

    def __init__(self, prcg):
        self.prcg = prcg

    @property
    def depends(self):
        yield self.dep_class(self.prcg)

    @property
    def fn(self):
        return ("{outdir}/{prcg:dir}/{code}/{prc:gen}.{fext}"
                .format(outdir=config.outdir,
                        code=self.code,
                        prc=self.prcg,
                        fext=self.fext)
                )

    @property
    def log_fn(self):
        return ("{outdir}/{prcg:dir}/{code}/{prcg:gen}.{fext}"
                .format(outdir=config.outdir,
                        code=self.code,
                        prcg=self.prcg,
                        fext='log')
                )

    @property
    def is_done(self):
        return os.path.exists(self.fn)

    def do(self, task):
        in_fn = task.fn  # isn't task going to be a list? why does this not throw exception?
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
        return "<{} {}>".format(self.__class__.__name__, self.prcg)

    def do_file(self, infn, outfn, logfn=None):
        raise NotImplementedError


class ProjRunClone(Task):
    dep_class = PRCGTask

    def __init__(self, project, run, clone, indir):
        self.project = project
        self.run = run
        self.clone = clone
        self.indir = indir
        self._depends = None

    def _configure(self, prcg):
        return prcg

    def get_gens(self, prc_dir):
        for fn in glob.iglob("{prc_dir}/*".format(prc_dir=prc_dir)):
            yield 0, fn

    def _get_prcgs(self):
        for gen, rawfn in self.get_gens(self.indir):
            yield self._configure(
                PRCG(self.project, self.run, self.clone, gen, rawfn)
            )

    @property
    def depends(self):
        if self._depends is None:
            self._depends = list(
                self.dep_class(prcg)
                for prcg in self._get_prcgs()
            )
        yield from self._depends

    def do(self, tasks):
        pass

    def is_done(self):
        pass


class Project(Dummy, Task):
    dep_class = ProjRunClone

    def __init__(self, project):
        self.project = project
        self.indir = ("{indir}/{project}"
                      .format(indir=config.indir, project=project))
        if not os.path.exists(self.indir):
            raise ValueError("Project in directory doesn't exist. "
                             "Looking for {}".format(self.indir))

        self._depends = None

    def get_run_clones(self, indir):
        for fn in glob.iglob("{indir}/*".format(indir=indir)):
            yield 0, 0, fn

    @property
    def depends(self):
        if self._depends is None:
            self._depends = list(
                self.dep_class(self.project, run, clone, prc_dir)
                for run, clone, prc_dir in self.get_run_clones(self.indir))
        yield from self._depends


class FahProject(Project):
    prc_glob = "{indir}/RUN*/CLONE*/"
    prc_re = r"{indir}/RUN(\d+)/CLONE(\d+)/"

    def get_run_clones(self, indir):
        for fn in glob.iglob(self.prc_glob.format(indir=indir)):
            ma = re.match(self.prc_re.format(indir=indir), fn)
            yield int(ma.group(1)), int(ma.group(2)), fn


class FahProjRunClone(ProjRunClone):
    gen_re = re.compile("")
    gen_glob = ""

    def get_gens(self, prc_dir):
        for fn in (glob.iglob(self.gen_glob.format(prc_dir=prc_dir))):
            yield int(self.gen_re.search(fn).group(1)), fn


class StructPerRun:
    """Mix this in to ProjRunClone task to add struct information

    This looks for a json file in the input directory of the form
    {project}-structs.json
    which is keyed by the run and has values "struct" and "fext"
    """

    def _configure(self, prcg):
        prcg = super()._configure(prcg)
        if not hasattr(self, 'structs'):
            with open("{indir}/{prcg.project}-structs.json"
                              .format(indir=config.indir, prcg=prcg)) as f:
                self.structs = json.load(f)

        prcg.meta['struct'] = self.structs[str(prcg.run)]['struct']
        prcg.meta['top_fext'] = self.structs[str(prcg.run)]['fext']
        prcg.meta['top_dir'] = ("{indir}/{prcg.project}-tops/"
                                .format(indir=config.indir, prcg=prcg))
        prcg.meta['top_fn'] = ("{prcg.meta[top_dir]}/"
                               "{prcg.meta[struct]}.{prcg.meta[top_fext]}"
                               .format(indir=config.indir, prcg=prcg))
        return prcg


class ProjRunClonex21(FahProjRunClone):
    gen_re = re.compile(r"results-(\d\d\d)/")
    gen_glob = "{prc_dir}/results-???/positions.xtc"


class ProjRunClonexA4(FahProjRunClone):
    gen_re = re.compile(r"frame(\d+).xtc")
    gen_glob = "{prc_dir}/frame*.xtc"

    def _configure(self, prcg):
        prcg = super()._configure(prcg)
        prcg.meta['needs_trjconv'] = True
        prcg.meta['has_overlapping_frames'] = True
        prcg.meta['tpr_fn'] = ("{indir}/frame0.tpr"
                               .format(indir=os.path.dirname(prcg.in_fn)))
        return prcg


class ProjectBluewaters(Project):
    def _configure(self, prcg):
        prcg.meta['needs_trjconv'] = True
        prcg.meta['tpr_fn'] = ("{indir}/topol.tpr"
                               .format(indir=os.path.dirname(prcg.in_fn)))
        return prcg

    def get_run_clones(self, indir):
        # TODO
        for fn in glob.iglob("*/"):
            yield 0, 0, fn

    def get_gens(self, prc_dir):
        # TODO
        for fn in glob.iglob("{prc_dir}/*.xtc".format(prc_dir=prc_dir)):
            yield 0, fn
