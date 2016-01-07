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
        yield PRC(project, int(ma.group(1)), int(ma.group(2))), fn


def get_incomplete_prcs(project, projtype, indir):
    if projtype in ['xa4', 'x21']:
        yield from get_prcs_fah(indir, project)


def get_gens_xa4(prc_indir):
    gen_re = re.compile(r"frame(\d+).xtc")
    for fn in glob.iglob("{prc_indir}/frame*.xtc".format(prc_indir=prc_indir)):
        yield int(gen_re.search(fn).group(1)), fn


def get_gens_x21(prc_indir):
    gen_re = re.compile(r"results-(\d\d\d)/")
    for fn in (glob.iglob("{prc_indir}/results-???/positions.xtc"
                                  .format(prc_indir=prc_indir))):
        yield int(gen_re.search(fn).group(1)), fn


def get_gens(prc, projtype, prc_indir):
    if projtype == 'xa4':
        yield from get_gens_xa4(prc_indir)
    elif projtype == 'x21':
        yield from get_gens_x21(prc_indir)


def get_prcs(project, projtype, indir):
    for prc1, prc_indir in get_incomplete_prcs(project, projtype, indir):
        for gen, rawfn in get_gens(prc1, projtype, prc_indir):
            prc = PRC(prc1.project, prc1.run, prc1.clone, gen, rawfn)
            if projtype in ['xa4', 'bw']:
                prc.flags.add("needs_trjconv")
            yield prc
