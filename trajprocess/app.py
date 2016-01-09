from contextlib import contextmanager
from datetime import datetime
import traceback
import sys


class MockLoadBalancedView:
    def map_async(self, func, args):
        for a in args:
            func(a)

    @contextmanager
    def temp_flags(self, after):
        yield


def dump(task, exc):
    fn = "./Exception-{}.log".format(datetime.now().isoformat())
    with open(fn, 'a') as f:
        f.write("\n".join([
            "An exception occured",
            "--------------------",
            "{}".format(task.__class__),
            "{}".format(task.prc) if hasattr(task, 'prc') else "unknown prc",
            "",
            str(task.__dict__),
            "",
            "",
            str(exc),
            "",
            "",
        ]))
        _, _, tb = sys.exc_info()
        traceback.print_tb(tb, file=f)


def _run_function(task):
    def _inner_run_function(depend):
        try:
            task.do(depend)
            return str(task)
        except Exception as exc:
            dump(task, exc)

    return _inner_run_function


def _execute(task, lbv):
    ars = []
    for dep in task.depends:
        if dep.is_ephemeral and task.is_done:
            pass
        else:
            ar = _execute(dep, lbv)
            if ar is not None:
                ars.append(ar)

    if not task.is_done:
        with lbv.temp_flags(after=ars):
            return lbv.map_async(_run_function(task), list(task.depends))


def execute_task(task, lbv=None):
    if lbv is None:
        lbv = MockLoadBalancedView()

    _execute(task, lbv)
