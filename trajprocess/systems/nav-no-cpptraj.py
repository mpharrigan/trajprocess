from .. import tasks
from ..process import run_trjconv, convert_to_nc


class Trjconv(tasks.PRCTask):
    """Run gmx trjconv on system simulated with gromacs mdrun

    mdrun doesn't reimage *at all*. This will exploit the tpr to image
    things as well as possible.
    """
    code = 'cnv1'
    fext = 'xtc'
    dep_class = tasks.RawXTC
    needs_log = True
    is_ephemeral = True

    def do_file(self, infn, outfn, logfn=None):
        assert self.prc.project == 'p9752', "stride 4 makes no sense otherwise."
        run_trjconv(
                infn, outfn,
                log_fn=logfn,
                top_fn=self.prc.meta['tpr_fn'],
                stride=4,
        )


class ConvertToNC(tasks.PRCTask):
    """Convert trajectories to Amber NetCDF

    cpptraj can't read gromacs xtc files.

    Note: This task will conditionally depend on Trjconv depending on
    prc metadata.
    """
    code = 'cnv2'
    is_ephemeral = True

    @property
    def depends(self):
        if "needs_trjconv" in self.prc.meta:
            yield Trjconv(self.prc)
        else:
            yield tasks.RawXTC(self.prc)

    def do_file(self, infn, outfn, logfn=None):
        overlap = 'has_overlapping_frames' in self.prc.meta
        convert_to_nc(infn, outfn, has_overlapping_frames=overlap)


class Clean(tasks.Clean):
    ephemeral_task_classes = [
        Trjconv,
    ]
    dep_class = ConvertToNC
    delete_empty_dirs = True
    delete_logs = True


class Projectx21(tasks.StructPerRun, tasks.Projectx21):
    dep_class = Clean


class ProjectxA4(tasks.StructPerRun, tasks.ProjectxA4):
    dep_class = Clean


class NaV(tasks.Dummy, tasks.Task):
    depends = [
        Projectx21("p9704"),
        ProjectxA4("p9752"),
    ]
