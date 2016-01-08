

class MockLoadBalancedView:
    def apply_async(self, func, args):
        for a in args:
            func(a)


def _execute(task, lbv):
    for dep in task.depends:
        _execute(dep, lbv)

    if not task.is_done:
        lbv.apply_async(task.do, task.depends)


def execute_task(task, lbv=None):
    if lbv is None:
        lbv = MockLoadBalancedView()

    _execute(task, lbv)
