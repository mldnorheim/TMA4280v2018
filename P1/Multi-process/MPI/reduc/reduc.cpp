#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <mpi.h>
#include <cstdio>
#include <string.h>
#include <cstdlib>
#include <cmath>
#include <numeric>
#include <functional>

using namespace std;

int main(int argc, char *argv[])
{
    int apps = atoi(argv[1]);
    int type = atoi(argv[2]);

    MPI_Init(&argc, &argv);



    // --------------------- ZETA ---------------------------


    if(apps == 0)
    {
        int nprocs, myrank;
        const double pi = 3.141592653589793;


        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        double time_start;
        int n;
        int N = log2(nprocs);  // nprocs = 2^N

        if (myrank == 0) {
            cout << "Please enter a value for the number of approx. terms, n: " << endl;
            cin >> n;
            time_start = MPI_Wtime();
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&n, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


        // MPI_Allreduce

        if(type==0)
        {
            double vR = 0.0;
            for(int j=myrank+1; j<=n; j+=nprocs){
                vR += pow(j,-2.0);
            }
            double SR;
            MPI_Allreduce(&vR, &SR, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            double piR = sqrt(6*SR);
            double error;

            if (myrank==0) {
                double duration = MPI_Wtime() - time_start;
                double error = fabs(piR-pi);
                printf(" piR=%e \n errorR=%e \n durationR=%e \n", piR, error, duration);
            }
        }

        // "recursive-duoubling sum" (method from Question 1)

        if(type == 1)
        {
            double sigmaR = 0.0;
            for(int j=myrank+1; j<=n; j+=nprocs){
                double vj = pow(j,-2.0);
                sigmaR += vj;
            }

            for (int d=0; d<=(N-1); d++) {
                double sigmaRq;
                MPI_Send(&sigmaR, 1, MPI_DOUBLE, myrank^(int)pow(2,d), 0, MPI_COMM_WORLD);
                MPI_Recv(&sigmaRq, 1, MPI_DOUBLE, myrank^(int)pow(2,d), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                sigmaR += sigmaRq;
            }

            double piR = sqrt(6.0*sigmaR);

            double error;

            if (myrank==0) {
                double duration = MPI_Wtime() - time_start;
                double error = fabs(piR-pi);
                printf(" piR=%e \n errorR=%e \n durationR=%e \n", piR, error, duration);
            }
        }

    }



    // --------------------- MACH ---------------------------


    if(apps == 1)
    {
        int nprocs, myrank;
        const double pi = 3.141592653589793;

        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        double time_start;
        int n;
        int N = log2(nprocs);  // nprocs = 2^N

        if (myrank == 0) {
            cout << "Please enter a value for the number of approx. terms, n: " << endl;
            cin >> n;
            time_start = MPI_Wtime();
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(&n, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


        // MPI_Allreduce

        if(type == 0)
        {
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
            MPI_Allreduce(&vM5, &SM5, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            MPI_Allreduce(&vM239, &SM239, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            double piM = 4.0*(4.0*SM5 - SM239);
            double error;
            if (myrank==0) {
                double duration = MPI_Wtime() - time_start;
                double error = fabs(piM-pi);
                printf(" piM=%e \n errorM=%e \n durationM=%e \n", piM, error, duration);
            }

        }


        // "recursive-duoubling sum"

        if(type == 1)
        {

            double sigmaM5 = 0.0;
            double sigmaM239 = 0.0;
            for(int j=1; j<=n; j++){
                // x=1/5
                double vjx5 = pow(-1.0,j-1) * pow(1.0/5.0,2.0*j-1) / (2.0*j-1.0);

                sigmaM5 += vjx5;

                // x=1/239
                double vjx239 = pow(-1.0,j-1) * pow(1.0/239.0,2.0*j-1) / (2.0*j-1.0);

                sigmaM239 += vjx239;
            }

            for (int d=0; d<=(N-1); d++) {

                double sigmaM5q;
                double sigmaM239q;
                MPI_Send(&sigmaM5, 1, MPI_DOUBLE, myrank^(int)pow(2,d), 0, MPI_COMM_WORLD);
                MPI_Recv(&sigmaM5q, 1, MPI_DOUBLE, myrank^(int)pow(2,d), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                sigmaM5 += sigmaM5q;
                MPI_Send(&sigmaM239, 1, MPI_DOUBLE, myrank^(int)pow(2,d), 0, MPI_COMM_WORLD);
                MPI_Recv(&sigmaM239q, 1, MPI_DOUBLE, myrank^(int)pow(2,d), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                sigmaM239 += sigmaM239q;

            }

            double piM = 4.0*(4.0*sigmaM5 - sigmaM239);
            double error;
            if (myrank==0) {
                double duration = MPI_Wtime() - time_start;
                double error = fabs(piM-pi);
                printf(" piM=%e \n errorM=%e \n durationM=%e \n", piM, error, duration);
            }

        }

    }






    MPI_Finalize();
    return 0;
}
