# Diffusion-Monte-Carlo
Diffusion Monte Carlo(DMC) Calculations on simple harmonic oscillator(SHO).  
This program is written in Fortran90 to calculate the ground state energy of SHO. For N dimensional case, the ground state energy should be 0.5N, with mass = omega = 1.0.
Actually I want to follow the algorithmic form given by J. M. Thijssen in his book "Computational Physics(Second Edition)". However, the updating of trial energy and the number of birth might have some problems, which result in the explosion of walkers. Therefore, I change the code of these part basing on other DMC code to make it works.

To use this code, just simply type as follows:
 ifort SHO-DMC.f90
 ./a.out
The results will be printed on the screen and written in file 'fort.77'.
