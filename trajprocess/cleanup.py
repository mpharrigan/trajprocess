from pathlib import Path
import json
import os
import logging

from .process import config

log = logging.getLogger(__name__)


def try_remove(info, step, i):
    try:
        os.remove(info[step]['gens'][i])
    except IndexError:
        log.debug("Index error {project}-{run}-{clone}"
                  .format(**info['meta']))
    except FileNotFoundError:
        log.debug("File does not exist {project}-{run}-{clone}"
                  .format(**info['meta']))


def clean(info):
    try:
        n_done = len(info['ctr']['gens'])
    except KeyError:
        n_done = 0
    log.info("Found {n_done} done in {project}-{run}-{clone}"
             .format(n_done=n_done, **info['meta']))
    for i in range(n_done):
        try_remove(info, 'cnv1', i)
        try_remove(info, 'cnv2', i)
        try_remove(info, 'stp', i)


def cleanup():
    home = Path(config.prefix)
    info_ps = home.glob("**/info.json")
    for info_p in info_ps:
        with info_p.open('r') as f:
            info = json.load(f)
            clean(info)
