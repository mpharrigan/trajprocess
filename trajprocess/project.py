import re
import glob

from .prc import PRC

SUPPORTED_TYPES = ['xa4', 'x21', 'bw']


def parse_project(project):
    ma = re.match(r"([a-z]+)([0-9]+)", project)
    if ma is None:
        raise ValueError("Invalid project name. "
                         "Must be of the form xx12345")
    return ma.group(1), int(ma.group(2))


def parse_projtype(type):
    if type not in SUPPORTED_TYPES:
        raise ValueError("Invalid project type: {}. "
                         "Must be one of {}"
                         .format(type, SUPPORTED_TYPES))
    return type


def get_prcs_fah(indir, project):
    for fn in glob.iglob("{indir}/RUN*/CLONE*/".format(indir=indir)):
        ma = re.match(r"{indir}/RUN(\d+)/CLONE(\d+)/"
                      .format(indir=indir), fn)
        yield PRC(project, int(ma.group(1)), int(ma.group(2)), fn)


def get_prcs(project, projtype, indir):
    if projtype in ['xa4', 'x21']:
        yield from get_prcs_fah(indir, project)


def get_gens(prc, projtype):
    if projtype == 'xa4':
        gen_re = re.compile(r"frame(\d+).xtc")
        yield from (int(gen_re.search(fn).group(1))
                    for fn in glob.iglob("{prc:indir}/frame*.xtc"
                                         .format(prc=prc)))
    elif projtype == 'x21':
        gen_re = re.compile(r"results-(\d\d\d)/")
        yield from (int(gen_re.search(fn).group(1))
                    for fn in glob.iglob("{prc:indir}/results-???/"
                                         .format(prc=prc)))
