gm4poisson
==========

This is a geometric multigrid solver for 1d, 2d, 3d Possion problems. 

The solver is free of matrix.

It could be paralleled by 2 or 4 color G-S iteration

========================================================================
========================================================================

There are several problems in the code.

1. I am not very familiar with writing makefile, where the makefile in the lib might not work.
2. I am now working on adding Openmp into the solver, which requires very careful plan to keep the balance between the threads. So I will continue focusing on it.
3. Profiling. To reach a better performance, I will look up the hotspots and optimized in the future.
4. Adding comments and generating doxygen.
5. Preconditioned CG are not okay, convergent rate is not as high as what is proved to be.

Finally, about the algorithm itself, there are two things need to be done:
1. More smoother.
2. Much more important. Now, the nx, ny, nz must be equal to 2^k. However, in the real problems, the real grids are always not satisfied with this condition. So I have a basic edition where nx, ny, nz could be any positive integer, and I will add them soon.

