#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <mpi.h>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <cmath>
#include <fstream>
#include <numeric>
#include <functional>
//#include "verification.h"


using namespace std;

int main(int argc, char **argv)      // MACHIN:
{

    ofstream out_data1("errmach1.txt");
    ofstream out_data2("durmach1.txt");


    int nprocs, myrank;
    const double pi = 3.141592653589793;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    for (int k=1; k<=24; k++)
    {


        double time_start;
        int n;
        if (myrank == 0) {

            n = pow(2,k);
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
        double SM5;               //= accumulate(vM5.begin(), vM5.end(), 0.0);
        double SM239;             //= accumulate(vM239.begin(), vM239.end(), 0.0);

        MPI_Reduce(&vM5, &SM5, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&vM239, &SM239, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        double piM = 4.0*(4.0*SM5 - SM239);
        double error;
        if (myrank==0) {
            double duration = MPI_Wtime() - time_start;
            double error = fabs(piM-pi);
            out_data1 << k << " " << scientific << setprecision(18) << log(error) << endl;
            out_data2 << k << " " << scientific << setprecision(18) << duration << endl;
        }


    }
    MPI_Finalize();
    return 0;

}

