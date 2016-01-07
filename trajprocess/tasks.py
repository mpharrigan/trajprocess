import os
import subprocess
import re
import glob


class Task:
    @property
    def is_done(self):
        return False

    def do(self, tasks):
        return


class config:
    indir = 'trajprocess.in.4'
    outdir = 'trajprocess.out.4'


class PRC:
    def __init__(self, project, run, clone, indir):
        self.project = project
        self.run = run
        self.clone = clone
        self.indir = indir

    @property
    def as_tuple(self):
        return self.project, self.run, self.clone

    def __format__(self, format_spec):
        if format_spec == 'dir':
            return "/".join(str(s) for s in self.as_tuple)
        elif format_spec == 'indir':
            return self.indir
        else:
            return "-".join(str(s) for s in self.as_tuple)


class Strip(Task):
    def __init__(self, prc, gen):
        self.prc = prc
        self.gen = gen

    pass


class Center(Task):
    def __init__(self, prc, gen):
        self.prc = prc
        self.gen = gen

    @property
    def depends(self):
        yield from (Strip(self.prc, self.gen),)

    @property
    def is_done(self):
        return os.path.exists("{outdir}/{prc:dir}/ctr/{gen}.nc"
                              .format(outdir=config.outdir,
                                      prc=self.prc,
                                      gen=self.gen))

    def do(self, strip):
        subprocess.check_call(
                ['touch',
                 "{outdir}/{prc:dir}/ctr/{gen}.nc"
                     .format(outdir=config.outdir,
                             prc=self.prc,
                             gen=self.gen)
                 ],
        )


class Project(Task):
    supported_types = ['xa4', 'x21', 'bw']

    def __init__(self, project, type):
        self.project = project

        ma = re.match(r"([a-z]+)([0-9]+)", project)
        if ma is None:
            raise ValueError("Invalid project name. "
                             "Must be of the form xx12345")
        self.projcode = ma.group(1)
        self.projnum = int(ma.group(2))

        if type not in self.supported_types:
            raise ValueError("Invalid project type: {}. "
                             "Must be one of {}"
                             .format(type, self.supported_types))
        self.projtype = type

        self.indir = "{indir}/{project}".format(indir=config.indir,
                                                project=project)

        if not os.path.exists(self.indir):
            raise ValueError("Project in directory doesn't exist. "
                             "Looking for {}".format(self.indir))

        self._depends = None

    def _get_prcs_fah(self):
        def prc_from_fn(fn):
            ma = re.match(r"{indir}/RUN(\d+)/CLONE(\d+)/"
                          .format(indir=self.indir), fn)
            return PRC(self.project, int(ma.group(1)), int(ma.group(2)), fn)

        yield from (prc_from_fn(fn)
                    for fn in glob.iglob("{indir}/RUN*/CLONE*/"
                                         .format(indir=self.indir)))

    def _get_prcs(self):
        if self.projtype in ['xa4', 'x21']:
            yield from self._get_prcs_fah()

    def _get_gens(self, prc):
        if self.projtype == 'xa4':
            gen_re = re.compile(r"frame(\d+).xtc")
            yield from (int(gen_re.search(fn).match(1))
                        for fn in glob.iglob("{prc:indir}/frame*.xtc"
                                             .format(prc=prc)))
        elif self.projtype == 'x21':
            gen_re = re.compile(r"results-(\d\d\d)/")
            yield from (int(gen_re.search(fn).match(1))
                        for fn in glob.iglob("{prc:indir}/results-???/"
                                             .format(prc=prc)))

    def _get_depends(self):
        for prc in self._get_prcs():
            for gen in self._get_gens(prc):
                yield Strip(prc, gen)

    @property
    def depends(self):
        if self._depends is None:
            self._depends = list(self._get_depends())

        yield from self._depends


class MockLoadBalancedView:
    def apply_async(self, func, args):
        for a in args:
            func(a)


def execute_task(task, lbv=None):
    if lbv is None:
        lbv = MockLoadBalancedView()

    all_dep_done = True
    for dependency in task.depends:
        all_dep_done = all_dep_done and dependency.is_done
        if not dependency.is_done:
            execute_task(dependency, lbv)

    if all_dep_done:
        return lbv.apply_async(task.do, task.depends)


class NaV(Task):
    @property
    def depends(self):
        yield from [
            Project("p9704", 'x21'),
            Project("p9752", 'xa4'),
        ]


def main():
    execute_task(NaV())
