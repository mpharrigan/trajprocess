from .. import tasks
from ..config import config
from ..postprocess import call_cpptraj_stp, call_cpptraj_ctr
from ..process import run_trjconv, convert_to_nc


class Trjconv(tasks.PRCGTask):
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
        assert self.prcg.project == 'p9752', "stride 4 makes no sense otherwise"
        run_trjconv(
                infn, outfn,
                log_fn=logfn,
                top_fn=self.prcg.meta['tpr_fn'],
                stride=4,
        )


class ConvertToNC(tasks.PRCGTask):
    """Convert trajectories to Amber NetCDF

    cpptraj can't read gromacs xtc files.

    Note: This task will conditionally depend on Trjconv depending on
    prcg metadata.
    """
    code = 'cnv2'
    is_ephemeral = True

    @property
    def depends(self):
        if "needs_trjconv" in self.prcg.meta:
            yield Trjconv(self.prcg)
        else:
            yield tasks.RawXTC(self.prcg)

    def do_file(self, infn, outfn, logfn=None):
        overlap = 'has_overlapping_frames' in self.prcg.meta
        convert_to_nc(infn, outfn, has_overlapping_frames=overlap)


class Strip(tasks.PRCGTask):
    """Use cpptraj to strip all but closest x."""
    code = 'stp'
    dep_class = ConvertToNC
    needs_log = True
    is_ephemeral = True

    def do_file(self, infn, outfn, logfn=None):
        call_cpptraj_stp(
                infn, outfn, logfn,
                removes=[":WAT", ":MY", "@Na+", "@Cl-"],
                num_to_keeps=[10000, 100, 20, 20],
                prmtopdir="{indir}/p9704-tops".format(indir=config.indir),
                outtopdir="{outdir}/prmtops".format(outdir=config.outdir),
                struct=self.prcg.meta['struct'],
        )


class Center(tasks.PRCGTask):
    """Use cpptraj to center and image."""
    code = 'ctr'
    dep_class = Strip
    needs_log = True

    def do_file(self, infn, outfn, logfn=None):
        call_cpptraj_ctr(
                infn, outfn, logfn,
                stptopdir="{outdir}/prmtops".format(outdir=config.outdir),
                struct=self.prcg.meta['struct'],
        )


class Clean(tasks.Clean):
    ephemeral_task_classes = [
        Strip,
        ConvertToNC,
        Trjconv
    ]
    dep_class = Center
    delete_logs = True


class PRCx21(tasks.StructPerRun, tasks.ProjRunClonex21):
    dep_class = Clean


class PRCxA4(tasks.StructPerRun, tasks.ProjRunClonexA4):
    dep_class = Clean


class Projectx21(tasks.FahProject):
    dep_class = PRCx21


class ProjectxA4(tasks.FahProject):
    dep_class = PRCxA4


class NaV(tasks.Dummy, tasks.Task):
    depends = [
        Projectx21("p9704"),
        ProjectxA4("p9752"),
    ]
