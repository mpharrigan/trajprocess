class PRCG:
    def __init__(self, project, run, clone, gen, in_fn=None):
        self.project = project
        self.run = run
        self.clone = clone
        self.gen = gen
        self.in_fn = in_fn
        self.meta = dict()

    @property
    def as_tuple(self):
        return self.project, self.run, self.clone

    def __format__(self, format_spec):
        if format_spec == 'dir':
            return "/".join(str(s) for s in self.as_tuple)
        if format_spec == 'gen':
            return str(self.gen)
        elif format_spec == 'raw':
            return self.in_fn
        else:
            return "-".join(str(s) for s in list(self.as_tuple) + [self.gen])
