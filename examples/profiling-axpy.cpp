#include "ptimer.hpp"
#include "mpi.h"

xns::ptimer chrono;

void axpy(int n, double a, double* x, double* y){
    chrono.start("axpy");
    /* manually count flops and memory accesses */
    chrono.add_flops(2*n);
    chrono.add_bytes(2*8*n);
    for(int i=0; i<n; ++i){
        y[i] = a*x[i] + y[i];
    }
    chrono.pause("axpy");
}

int main(int argc, char **argv){
    MPI_Init(&argc, &argv);

    const int N = 1e8;

    chrono.start("setup");
    chrono.start("allocation");
    double* x = new double[N];
    double* y = new double[N];
    chrono.pause("allocation");

    chrono.start("initialization");
    for(int i=0; i<N; ++i){
        x[i] = (double)i;
    }
    chrono.pause("initialization");
    chrono.pause("setup");

    chrono.start("compute");
    axpy(N, 1.2345, x, y);
    axpy(N, 0.1234, y, x);
    axpy(N, 0.0123, x, y);
    axpy(N, 0.0012, y, x);
    chrono.pause("compute");
    chrono.print();
    
    delete[] x;
    delete[] y;

    MPI_Finalize();

    return 0;
}
