# ptimer
Parallel timer for high-performance computing.

To run the example, use the following command:

FORTRAN
mpif90 ptimer.f90 ./examples/profiling-axpy.f90 -o ./examples/profiling-axpy.out; ./examples/profiling-axpy.out

C++
mpic++ -std=c++11 -I. ./examples/profiling-axpy.cpp -o ./examples/profiling-axpy.out; ./examples/profiling-axpy.out
