# `trajprocess`

Take MD trajectories from a variety of sources and normalize them.

Each trajectory can be in multiple chunks ("gens" in FAH vocabulary). The
processing on each gen will not be duplicated. The files will never be concatenated.
Use mdtraj's support for loading multiple files per trajectory when doing your
analysis: `md.load(["fn1.nc", "fn2.nc"], top="whatever.prmtop")`

An extensive cache of metadata is stored in one `info.json` file per trajectory.
Your analysis should start by grepping all the json files, saving the order
in which they were grepped, and loading trajectories based on the contents of
the json file: `md.load(info['ctr']['gens'], top=info['stp']['outtop'])`

This package was developed to support combining the output of heterogeneous
MD engines (openmm, gromacs, and amber), heterogeneous hpc environments
(supercomputers, fah), and slightly heterogeneous topologies
(different numbers of waters, lipids, ions). Output trajectories should be
immediately ready for any quantitative or visual analysis.

Trajectories are processed in stages. Each stage codename serves as a key
in the metadata `info.json`:

 - `cnv1`: optionally convert trajectories using `gmx trjconv -pbc whole`
 - `cnv2`: convert to netcdf using `mdtraj` for use in `cpptraj`
 - `stp`: strip all but closest N waters, lipids, ions using `cpptraj`
 - `ctr`: autoimage and center using `cpptraj`

Other metadata includes:

 - `raw`: input
 - `meta`: project, run, clone information
 - `path`: relevant filesystem paths

## Requirements

 - `mdtraj` >= 1.5 (fix for netcdf writing)
 - `gmx trjconv` from Gromacs
 - `cpptraj` from AmberTools


## Customization

Currently, you have to edit the source code to customize settings
for new projects or project types.
