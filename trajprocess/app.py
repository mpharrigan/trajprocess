class config:
    indir = 'trajprocess.in.4'
    outdir = 'trajprocess.out.4'


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
