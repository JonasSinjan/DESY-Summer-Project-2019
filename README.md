# DESY

Code files for Summer Student 2019 at DESY Zeuthen Campus. //
Working in the Astroparticle Physics Group led by Prof. Huirong Yan. //
Supervisor: Dr. Kirit Makwana //

Runtimes:

-Currently the 3d displacement 128 with Phi0 real initialisation requires 43 hours on a single thread.

-Currently the 2d displacement 512 with Phi0 real initialisation requires 40 minutes on a single thread.

-3d displacement now setup in a configuration which initialises Phi0 in fourier space: 128 takes 4 minutes.

-3d 512 requries approx 100GB RAM at least, needs to be made more memory efficient.

To Run:

1. Change the Makefile directory path for the fftw library that you must install locally, not present in all fortran compilers.

2. For all methods, change the number of grids variable: ngrids to the power of 2  of the resolution desired.

3. Change the output directory to the desired path.

4. Type 'make' in terminal to compile, then './run.x' in the method directory to execute the FORTRAN files.



