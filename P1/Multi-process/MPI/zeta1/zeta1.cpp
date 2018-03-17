#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <cmath>
#include <mpi.h>
#include <cstdio>
#include <string.h>
#include <cstdlib>
#include <math.h>
#include <numeric>
#include <functional>


using namespace std;

int main(int argc, char **argv)   // RIEMANN ZETA:
{

    int nprocs, myrank;
    const double pi = 3.141592653589793;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    double time_start;
    int n;

    if (myrank == 0) {
        cout << "Please enter a value for the number of approx. terms, n: " << endl;
        cin >> n;
        time_start = MPI_Wtime();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    double vR = 0.0;
    for(int j=myrank+1; j<=n; j+=nprocs){
        vR += pow(j,-2.0);
    }
    double SR;
    MPI_Reduce(&vR, &SR, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    double piR = sqrt(6*SR);
    double error;

    if (myrank==0) {
        double error = fabs(piR-pi);
        double duration = MPI_Wtime() - time_start;
        //cout << "MPI-Riemann approximates pi as: " << scientific << piR << endl;
        //cout << "The error is: " << scientific << errorR << endl;
        printf(" piR=%e \n errorR=%e \n durationR=%e \n", piR, error, duration);
    }
    MPI_Finalize();

    return 0;
}

