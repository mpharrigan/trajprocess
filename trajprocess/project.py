import re

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
