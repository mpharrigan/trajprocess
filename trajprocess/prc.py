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
