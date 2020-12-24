# Diffusion-Monte-Carlo-SHO
Diffusion Monte Carlo(DMC) Calculations of simple harmonic oscillator(SHO) system.  

This program is written in Fortran90 to calculate the ground state energy of SHO system. For N dimensional case, the ground state energy should be 0.5\*N, with mass = omega = 1.0.

To use this code, just simply type as follows:

 ifort SHO-DMC.f90
 
 ./a.out
 
The results will be printed on the screen and written in file 'SHO_DMC_results.log'.
