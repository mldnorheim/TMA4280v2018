#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <omp.h>
#include <cstdio>
#include <cstdlib>
#include <string.h>
#include <cmath>
#include <numeric>
#include <functional>
#include <fstream>

using namespace std;

int main(int argc, char **argv)
{

    if (argc < 2){
        printf("Lacks an argument: Remembered nthreads?\n");
        return 1;
    }

    int nthreads = atoi(argv[1]);
    const double pi = 3.141592653589793;

    ofstream out_data1("errmach2.txt");
    ofstream out_data2("durmach2.txt");

    for (int k=1; k<=24; k++)
    {
        int n = pow(2,k);
        double time_start = omp_get_wtime();

        double sumM5 = 0.0;
        double sumM239 = 0.0;

        #pragma omp parallel for num_threads(nthreads) reduction(+:sumM5,sumM239)
            for(int j=1; j<=n; j++){
                // x=1/5
                sumM5 += pow(-1.0,j-1) * pow(1.0/5.0,2.0*j-1) / (2.0*j-1.0);
                // x=1/239
                sumM239 += pow(-1.0,j-1) * pow(1.0/239.0,2.0*j-1) / (2.0*j-1.0);
            }

        double piM = 4.0*(4.0*sumM5 - sumM239);
        double duration = omp_get_wtime() - time_start;
        double error = fabs(piM-pi);

        out_data1 << k << " " << scientific << setprecision(18) << log(error) << endl;
        out_data2 << k << " " << scientific << setprecision(18) << duration << endl;
    }

    return 0;

}


