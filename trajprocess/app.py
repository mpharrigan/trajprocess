from contextlib import contextmanager
from datetime import datetime
import traceback
import sys
import logging

log = logging.getLogger(__name__)


class MockLoadBalancedView:
    def map_async(self, func, args):
        for a in args:
            try:
                func(a)
            except:
                pass

    @contextmanager
    def temp_flags(self, **kwargs):
        yield


def dump(task, exc):
    fn = "./Exception-{}.log".format(datetime.now().isoformat())
    with open(fn, 'a') as f:
        f.write("\n".join([
            "An exception occured",
            "--------------------",
            "{}".format(task),
            "{}".format(task.__class__),
            "{}".format(task.prcg) if hasattr(task, 'prcg') else "unknown prcg",
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


class _run_function:
    def __init__(self, task):
        self.task = task

    def __call__(self):
        try:
            self.task.do()
            return str(self.task)
        except Exception as exc:
            dump(self.task, exc)
            raise exc


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
        with lbv.temp_flags(after=ars, retries=10):
            log.info("Submitting {}".format(task))
            return lbv.apply_async(_run_function(task))


def execute_task(task, lbv=None):
    if lbv is None:
        lbv = MockLoadBalancedView()

    _execute(task, lbv)
