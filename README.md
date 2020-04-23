# Stokes Repository

Code to solve the Stokes equations using boundary integral equations.

## Features
* confined flows, uses Power-Miranda formulation and an FMM to rapidly compute sums
* periodic domains, uses spectral Ewald to rapidly compute infinte sums
* special quadrature for on-surface evaluation of single-layer potential as well as near-surface evaluation of single- and double-layer potentials

## Compilation instructions

In addition to this repository you will need the fast multiplication methods in http://github.com/lukasbystricky/FastTools2D. Compilation instructions for that repository are provided in its README, however the following sequence has been tested on Ubuntu 16.04 and Matlab 2017a:

	## clone repositories
	git clone http://github.com/lukasbystricky/FastTools2D.git
	git clone http://github.com/lukasbystricky/Stokes2D.git

	## compile periodic spectral Ewald
	cd FastTools2D/PeriodicEwald/src/build
	cmake ..
	make

	## compile FMM if desired
	## NB: on Mac, Matlab mex requires the intel compilers
	## The makefiles will need to be edited!
	cd FastTools2D/FMM/src/StokesDLP
	make
	cd ../StokesSLP
	make

	## compile special quadrature routines
	cd ../../../Stokes2D/src/build
	cmake ..
	make
	cd ../..


Start Matlab and run `initpaths.m` to initialize the proper paths. To test the spectral Ewald and special quadrature run `demo_periodic_pouseille_flow`. With 10 panels on each wall, you should acheive around 14 digits of accuracy everywhere in the domain.  

To test the FMM, run `demo_taylor_couette.m`. With 10 panels on the outer wall and 5 panels on the inner wall you should see 12 digits of accuracy everywhere.

## To do

* check OpenMP timings
* RCIP for corner treatment
* Navier-slip boundary conditions, adjoint double-layer potential
