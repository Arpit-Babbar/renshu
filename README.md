# arpit_practise

The latest C++ code completed on 13th January was fv2d_dirichlet which solves a coefficient/variable advection PDE with Dirichlet bc using Upwind/Lax-Wendroff. It is just a modification of the Periodic solver of same PDE in fv2d_var_coeff.

The makefile has its dependencies in a separate include directory. The reason to keep a separate directory is that both files use the same `vtk_anim, array2d, initial_conditions` files.

**To do**

Merge fv2d_var_coeff and fv2d_dirichlet into the same directory; make a run.cc for user to choose which one to run
