import multiprocessing
from contextlib import contextmanager


class Pool(multiprocessing.Pool):
    def map_async(self, func, iterable, chunksize=None, callback=None,
                  error_callback=None):
        super().map(func, iterable, chunksize=1)

    @contextmanager
    def temp_flags(self, **kwargs):
        yield
