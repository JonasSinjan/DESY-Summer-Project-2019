# DESY

Code files for Summer Research Student Project (2019) at DESY Zeuthen Campus. //
Working in the Astroparticle Physics Group led by Prof. Huirong Yan. //
Supervisor: Dr. Kirit Makwana //

Project: Creating synthetic MHD turbulence data in astrophysical plasmas

Motivation: Standard numerical simulations that solve the incompressible MHD equations directly take ~1 million CPU hours for a sizeable resolution - this project was to investigate if turbulence data could be created much more cheaply - turbulence data like this used for cosmic ray diffusion studies

Credit: Michael Vorster for developing the FORTRAN code for the 2d squares & displacment method as well as the initial structure function python code. 

My work was to get these running myself, testing them, extending 2d displacement to 3d, and running various tests to analyse problems and validate the synthetic data with python: structure function, contour plots, power spectra etc

Created two more FORTRAN methods: phi0init - used for testing the scalar field initialisation in fourier space
                                  , localB - to create local B fields with perturbations by adding waves, while ensuring div.B = 0

Runtimes:

-Currently the 3d displacement 128 with Phi0 real initialisation requires 43 hours on a single thread.

-Currently the 2d displacement 512 with Phi0 real initialisation requires 40 minutes on a single thread.

-3d displacement now setup in a configuration which initialises Phi0 in fourier space: 128 takes 4 minutes.

-3d 512 requries approx 100GB RAM at least, needs to be made more memory efficient - hence 3d_disp_mem method - only 60GB for 512

To Run:

1. Change the Makefile directory path for the fftw library that you must install locally, not present in all fortran compilers.

2. For all methods, change the number of grids variable: ngrids to the power of 2  of the resolution desired.

3. Change the output directory to the desired path.

4. Type 'make' in terminal to compile, then './run.x' in the method directory to execute the FORTRAN files.



