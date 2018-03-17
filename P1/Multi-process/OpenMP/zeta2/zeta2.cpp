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

    double sumR = 0.0;

    #pragma omp parallel for num_threads(nthreads) reduction(+:sumR)
        for(int j=1; j<=n; j++){
            sumR += 1.0/((double)j*(double)j);
        }

    double piR = sqrt(6*sumR);
    double error = fabs(piR-pi);
    double duration = omp_get_wtime() - time_start;

    printf(" piR=%e \n errorR=%e \n durationR=%e \n", piR, error, duration);

    return 0;


}




