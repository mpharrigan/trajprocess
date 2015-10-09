# `trajprocess`

Process trajectories in stages

 - `nfo`: generate metadata for each trajectory and save as json
 - `cat`: concatenate trajectories using `gmx trjcat`
 - `cnv`: optionally convert trajectories using `gmx trjconv -pbc whole`
 - `cnv` (part 2): convert to netcdf using `mdtraj`
 - `stp`: strip all but closest N waters, lipids, ions using `cpptraj`
 - `ctr`: autoimage and center using `cpptraj`
 
## Requirements

 - `mdtraj` with PR #939 (fix netcdf standard)
 - Gromacs >= 5.0.5 for bug #1705 (fix segfault during `trjcat` append)
