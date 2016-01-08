from ..tasks import Task, PRCTask, RawXTC, Projectx21 as _Projectx21, \
    ProjectxA4 as _ProjectxA4

from ..process import run_trjconv, convert_to_nc


class Trjconv(PRCTask):
    """Run gmx trjconv on system simulated with gromacs mdrun

    mdrun doesn't reimage *at all*. This will exploit the tpr to image
    things as well as possible.
    """
    code = 'cnv1'
    fext = 'xtc'
    dep_class = RawXTC
    needs_log = True

    def do_file(self, infn, outfn, logfn=None):
        assert self.prc.project == 'p9752', "stride 4 makes no sense otherwise."
        run_trjconv(
                infn, outfn,
                log_fn=logfn,
                top_fn=self.prc.meta['tpr_fn'],
                stride=4,
        )


class ConvertToNC(PRCTask):
    """Convert trajectories to Amber NetCDF

    cpptraj can't read gromacs xtc files.

    Note: This task will conditionally depend on Trjconv depending on
    prc metadata.
    """
    code = 'cnv2'

    @property
    def depends(self):
        if "needs_trjconv" in self.prc.meta:
            yield Trjconv(self.prc)
        else:
            yield RawXTC(self.prc)

    def do_file(self, infn, outfn, logfn=None):
        overlap = 'has_overlapping_frames' in self.prc.meta
        convert_to_nc(
                infn, outfn,
                has_overlapping_frames=overlap,
        )

        pass


class Strip(PRCTask):
    """Use cpptraj to strip all but closest x."""
    code = 'stp'
    dep_class = ConvertToNC

    def do_file(self, infn, outfn, logfn=None):
        pass


class Center(PRCTask):
    """Use cpptraj to center and image."""
    code = 'ctr'
    dep_class = Strip

    def do_file(self, infn, outfn, logfn=None):
        pass


class Projectx21(_Projectx21):
    dep_class = Center


class ProjectxA4(_ProjectxA4):
    dep_class = Center


class NaV(Task):
    depends = [
        Projectx21("p9704"),
        ProjectxA4("p9752"),
    ]
