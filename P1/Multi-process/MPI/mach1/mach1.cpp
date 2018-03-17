#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <mpi.h>
#include <cstdio>
#include <string.h>
#include <cmath>
#include <fstream>
#include <numeric>
#include <functional>
//#include "verification.h"


using namespace std;

int main(int argc, char **argv)      // MACHIN:
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

    double vM5 = 0.0;
    double vM239 = 0.0;
    for(int j=myrank+1; j<=n; j+=nprocs){
        // x=1/5
        vM5 += pow(-1.0,j-1) * pow(1.0/5.0,2.0*j-1) / (2.0*j-1.0);

        // x=1/239
        vM239 += pow(-1.0,j-1) * pow(1.0/239.0,2.0*j-1) / (2.0*j-1.0);
    }
    double SM5;
    double SM239;

    MPI_Reduce(&vM5, &SM5, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&vM239, &SM239, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    double piM = 4.0*(4.0*SM5 - SM239);
    double error;
    if (myrank==0) {
        double error = fabs(piM-pi);
        double duration = MPI_Wtime() - time_start;
        printf(" piM=%e \n errorM=%e \n durationM=%e \n", piM, error, duration);
        cout << "Machin approximates pi as: " << fixed << scientific << setprecision(15) << piM << endl;
        cout << "The error is: " << scientific << errorM << endl;
    }

    MPI_Finalize();

    return 0;

}

