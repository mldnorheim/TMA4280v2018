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

using namespace std;

int main(int argc, char **argv)
{

    if (argc < 3){
        printf("Need more arguments: Remembered nthreads and n?\n");
        return 1;
    }

    int nthreads = atoi(argv[1]);
    int n = atoi(argv[2]);

    const double pi = 3.141592653589793;
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

    printf(" piM=%e \n errorM=%e \n durationM=%e \n", piM, error, duration);
        //cout << "Machin approximates pi as: " << fixed << scientific << setprecision(15) << piM << endl;
        //cout << "The error is: " << scientific << errorM << endl;


    return 0;

}


