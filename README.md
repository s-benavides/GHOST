# GHOST

%%%%%
Addition notes by Santiago Benavides (santib@mit.edu)

All changes _only_ apply to solvers:
* HD
* ROTH
* MHD
* MHDB
* RMHDB

Changes:
* Added CFL condition.
* Added hyper- and hypo-dissipation.
* Modified 'initialfv.f90_random' to include the possibility of constant injection of energy (of amplitude f0) if 'rand' = 1. Added 'initialfv.f90_random_2', which is a different forcing function meant to be used for random forcing.
* Modified some of the global outputs. Including $|u(k_f)|^2$, $|J\times B|^2$ and hybrid helicity for RMHDB when the magnetic field and rotation are aligned.
* For RMHDB, added two new functions (with respective outputs): 
	(1) 'spec2D_yavg', which takes the 2D spectrum, but instead of integrating azimuthally, it integrates in the y-direction. Meant to be used with misaligned $\Omega$ and $B_0$ (with $B_0$ only x and z, not y component).
	(2) 'energy_arbdir', measures the total kinetic energy, split up into a parallel and perpendicular component, with respect to an arbitrary vector that is an input, e.g. $B_0$ or $\Omega".

%%%%%

GHOST (the Geophysical High-Order Suite for Turbulence) is an accurate and highly scalable pseudospectral code that solves a variety of PDEs often encountered in studies of turbulent flows. It is mainly developed by Pablo Mininni and Duane Rosenberg, with contributions from many users. The code uses a hybrid parallelization method combining MPI and OpenMP (GHOST also has support for GPUs using CUDA). This allows the code to run efficiently in laptops and small clusters, as well as to scale up to over 100,000 CPU cores in production runs in supercomputers. The hybrid parallelization method used in GHOST was recognized by two awards: an NCAR/CISL award, and the best paper award (technology track) at the TeraGrid 2010 conference.

The stable release found in the "master" branch supports cubic periodic domains with NxNxN grid points (in three-dimensions). For elongated domains with arbitrary linear sizes, and arbitrary resolutions in each spatial direction, check out the "pre-release-2" branch of GHOST.

GHOST can solve PDEs in periodic domains in two and three dimensions, to tackle many problems including:

* Compressible and incompressible hydrodynamic equations.
* Compressible and incompressible magnetohydrodynamic (MHD) equations.
* Hall-MHD equations and other kinetic plasma effects.
* Rotating flows.
* Rotating and stratified flows using the Boussinesq approximation.
* Passive scalars.
* Superfluids and Bose-Einstein condensates using the Gross-Pitaevskii Equation (GPE).
* Ginzburg-Landau equations.
* Integration of many types of particles.
* Several subgrid scale (SGS) models.

GHOST can integrate trajectories of Lagrangian particles, test particles (charged particles in a conducting fluid or a plasma), and inertial particles. It also includes two types of SGS models: one group of models derived from regularizations of the primitive equations of motion (Lagrangian Averaged, Clark, and Leray models), and another group of models based on computation of eddy viscosity and eddy noise using the Eddy-Damped Quasi-Normal Markovian approximation.

GHOST is mostly written in Fortran 90/95, with some C bindings, code in Fortran 2003 for the particles, and code in CUDA for GPUs. To build the main code, without particles and without GPU support, only a C and a Fortran 95 compiler are required.  More details on how to compile GHOST can be found in the directory 3D/src/README.txt. GHOST also includes examples of scripts in Python, Matlab, IDL, and Mathematica, to do post-analysis of the numerical simulations. Moreover, it can be used together with VAPOR (the Visualization and Analysis Platform for Ocean, Atmosphere, and Solar Researchers) to do still frame renderings and animations.  VAPOR provides an interactive environment for 3D visualization of GHOST outputs, allowing powerful analysis of turbulent flows in systems with 3D graphics cards (https://www.vapor.ucar.edu/).
