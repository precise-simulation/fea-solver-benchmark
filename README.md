Fortan, Julia, Octave, and Matlab FEM Benchmark and Comparison
==============================================================

Benchmark and comparison of Fortran, Julia, Octave, and Matlab for a
FEA Poisson problem solved on a unit square. The problem is
discretized with Q1 bilinear Lagrange finite elements.

The Fortran code is using a reference FEM implementation from which
the Julia code is a direct port. The Octave and Matlab code is derived
from the [FEATool](https://www.featool.com)
[Multiphysics](https://www.featool.com/multiphysics) production
code. The problem setup is identical and the codes are equivalent and
not tuned to perform better in the benchmark, as they still feature
all relevant code paths to allow different finite element shape
functions, variable coefficients, unstructured grids, etc. Thus the
results are relevant and comparable. The main difference is that
Octave, Matlab, and Julia use the default sparse linear solver
(currently the Umfpack direct solver), while the Fortran code uses a
more efficient iterative geometric multigrid solver.

The results of this benchmark are presented on the FEATool blog:
[https://www.featool.com/fortran-julia-and-matlab-fem-benchmark-comparison](https://www.featool.com/fem/2016/10/13/fortran-julia-and-matlab-fem-benchmark-comparison)

Installation and Running
------------------------

- Download and unzip or clone the repository.

- Set up and install [Octave](https://www.gnu.org/software/octave),
  Matlab [Julia](http://julialang.org), and a Fortran compiler
  (tested with gfortran and the Intel Fortran compiler).

- [Download and install the FEATool toolbox](https://www.featool.com/download)

- Edit the [testrun_param.txt](https://github.com/precise-simulation/fea-solver-benchmark/blob/master/testrun_param.txt#L1)
  file, which contains three parameters

        N0    - the grid resolution of the coarsest test grid (grid level)
                (the number of cells in x and y-directions)
        NLEV  - The number of grid levels to run tests for
        NRUNS - Number of test runs for each grid level,
                the timings are averaged for all runs

- On Windows edit the _OCTAVE_, _MATLAB_, and _JULIA_ paths in the
  [run_tests.bat](https://github.com/precise-simulation/fea-solver-benchmark/blob/master/run_tests.bat#L6)
  script and execute the script to automatically run the benchmarks
  and generate the output files. (Note that to run the Fortran code
  under Windows the _Windows Subsystem for Linux (WSL)_ is required)

- On other systems the **run_matlab.m**, **run_julia.jl**, and
  **run_fortran.sh** shell scripts can be run manually.

- The Octave/Matlab postprocessing script
  **src_matlab/process_results.m** file can be run manually
  to generate the results and output files.

- The [NNWORK](https://github.com/precise-simulation/fea-solver-benchmark/blob/master/src_fortran/src/featfem.f#L9)
  parameter in the main Fortran source file
  *src_fortran/src/featfem.f* controls the static memory allocation
  and might have to be increased and recompiled to run > 4 GB runs.


Results and Output
------------------

After running the tests the _output_ directory will contain results
files, _table.txt_ with tabulated times, and
[_results.html_](http://htmlpreview.github.io/?https://github.com/precise-simulation/fea-solver-benchmark/blob/master/output/results.html)
which shows the comparison graphs. Sample output is shown here below

![Assembly Timings](https://github.com/precise-simulation/fea-solver-benchmark/blob/main/output/fig_assembly.jpg)

![Solver Timings](https://github.com/precise-simulation/fea-solver-benchmark/blob/main/output/fig_solve.jpg)

![Total Timings](https://github.com/precise-simulation/fea-solver-benchmark/blob/main/output/fig_total.jpg)

Data for the complete results for all testruns are tabulated below

    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|
    |  Octave |           |           |           |           |           |           |           |           |           |
    |     1/h |    t_grid |     t_ptr |   t_asm_A |   t_asm_f |     t_bdr |  t_sparse |     t_tot |    t_spmv |   t_solve |
    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|
    |      32 | 2.37e-003 | 1.50e-003 | 1.77e-002 | 1.07e-002 | 1.06e-003 | 4.16e-004 | 3.37e-002 | 4.62e-005 | 5.78e-003 |
    |      64 | 2.76e-003 | 1.95e-003 | 1.99e-002 | 1.15e-002 | 2.76e-003 | 1.84e-003 | 4.07e-002 | 1.55e-004 | 2.28e-002 |
    |     128 | 6.28e-003 | 5.33e-003 | 4.39e-002 | 2.27e-002 | 1.13e-002 | 9.03e-003 | 9.86e-002 | 6.29e-004 | 1.04e-001 |
    |     256 | 1.74e-002 | 2.43e-002 | 1.28e-001 | 7.38e-002 | 4.25e-002 | 3.30e-002 | 3.19e-001 | 2.38e-003 | 5.10e-001 |
    |     512 | 7.38e-002 | 8.25e-002 | 4.92e-001 | 2.73e-001 | 1.57e-001 | 1.35e-001 | 1.21e+000 | 1.07e-002 | 2.66e+000 |
    |    1024 | 2.87e-001 | 3.28e-001 | 1.88e+000 | 1.07e+000 | 6.69e-001 | 5.33e-001 | 4.77e+000 | 4.33e-002 | 1.75e+001 |
    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|

    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|
    |  Matlab |           |           |           |           |           |           |           |           |           |
    |     1/h |    t_grid |     t_ptr |   t_asm_A |   t_asm_f |     t_bdr |  t_sparse |     t_tot |    t_spmv |   t_solve |
    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|
    |      32 | 1.74e-003 | 1.26e-003 | 2.50e-003 | 1.88e-003 | 2.22e-003 | 2.05e-003 | 1.16e-002 | 1.18e-005 | 4.05e-003 |
    |      64 | 7.09e-004 | 9.97e-004 | 6.51e-003 | 3.00e-003 | 6.13e-003 | 9.71e-003 | 2.71e-002 | 4.62e-005 | 2.17e-002 |
    |     128 | 1.58e-003 | 6.25e-003 | 1.58e-002 | 8.54e-003 | 2.28e-002 | 4.07e-002 | 9.57e-002 | 1.84e-004 | 9.43e-002 |
    |     256 | 5.06e-003 | 2.44e-002 | 5.49e-002 | 4.11e-002 | 9.02e-002 | 1.47e-001 | 3.63e-001 | 8.68e-004 | 4.48e-001 |
    |     512 | 3.06e-002 | 9.75e-002 | 2.05e-001 | 1.57e-001 | 3.74e-001 | 6.11e-001 | 1.47e+000 | 4.90e-003 | 2.54e+000 |
    |    1024 | 1.23e-001 | 3.86e-001 | 8.13e-001 | 6.32e-001 | 1.51e+000 | 2.51e+000 | 5.97e+000 | 1.99e-002 | 1.48e+001 |
    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|

    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|
    |   Julia |           |           |           |           |           |           |           |           |           |
    |     1/h |    t_grid |     t_ptr |   t_asm_A |   t_asm_f |     t_bdr |  t_sparse |     t_tot |    t_spmv |   t_solve |
    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|
    |      32 | 1.11e-004 | 3.77e-004 | 8.97e-004 | 2.61e-004 | 6.76e-005 | 2.30e-004 | 1.94e-003 | 2.15e-005 | 6.00e-003 |
    |      64 | 2.22e-004 | 1.67e-003 | 3.58e-003 | 1.00e-003 | 9.88e-005 | 1.20e-003 | 7.78e-003 | 8.86e-005 | 2.80e-002 |
    |     128 | 5.66e-004 | 6.76e-003 | 1.45e-002 | 3.92e-003 | 1.41e-004 | 5.96e-003 | 3.19e-002 | 3.65e-004 | 1.20e-001 |
    |     256 | 3.20e-003 | 2.71e-002 | 5.79e-002 | 1.58e-002 | 2.48e-004 | 2.39e-002 | 1.28e-001 | 1.55e-003 | 5.12e-001 |
    |     512 | 1.30e-002 | 1.10e-001 | 2.35e-001 | 6.66e-002 | 4.11e-004 | 1.77e-001 | 6.02e-001 | 7.77e-003 | 2.58e+000 |
    |    1024 | 5.35e-002 | 4.60e-001 | 9.57e-001 | 2.76e-001 | 2.22e-003 | 6.20e-001 | 2.37e+000 | 3.12e-002 | 1.34e+001 |
    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|

    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|
    | Fortran |           |           |           |           |           |           |           |           |           |
    |     1/h |    t_grid |     t_ptr |   t_asm_A |   t_asm_f |     t_bdr |  t_sparse |     t_tot |    t_spmv |   t_solve |
    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|
    |      32 | 0.00e+000 | 8.68e-004 | 8.68e-004 | 0.00e+000 | 0.00e+000 | 0.00e+000 | 1.74e-003 | 2.60e-005 | 4.34e-003 |
    |      64 | 0.00e+000 | 8.68e-004 | 7.81e-003 | 2.60e-003 | 8.68e-004 | 0.00e+000 | 1.22e-002 | 9.55e-005 | 8.68e-004 |
    |     128 | 4.34e-003 | 4.34e-003 | 1.82e-002 | 1.74e-003 | 0.00e+000 | 0.00e+000 | 2.86e-002 | 2.95e-004 | 1.13e-002 |
    |     256 | 3.47e-003 | 1.91e-002 | 6.94e-002 | 1.13e-002 | 1.74e-003 | 0.00e+000 | 1.05e-001 | 1.12e-003 | 4.08e-002 |
    |     512 | 3.30e-002 | 9.20e-002 | 2.26e-001 | 4.17e-002 | 1.74e-003 | 0.00e+000 | 3.94e-001 | 4.93e-003 | 1.77e-001 |
    |    1024 | 1.35e-001 | 3.86e-001 | 9.05e-001 | 1.67e-001 | 8.68e-004 | 0.00e+000 | 1.59e+000 | 2.35e-002 | 7.65e-001 |
    |---------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------|



License
-------

Copyright (C) 2013-2022 Precise Simulation Ltd.

Keywords: Finite Element, FEA, FEM, Fortran,
Julia, Octave, Matlab, Benchmark

This program is free software; you can redistribute it and/or modify
it under the terms of version 3 of the GNU Affero General Public
License (AGPLv3) as published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU Affero General Public
License along with this program. If not, see
[http://www.gnu.org/licenses](http://www.gnu.org/licenses).
