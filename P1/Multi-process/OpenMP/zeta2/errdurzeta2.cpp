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

    ofstream out_data1("errzeta2P8log.txt");
    ofstream out_data2("durzeta2P8.txt");

    int nthreads = atoi(argv[1]);
    const double pi = 3.141592653589793;

    for (int k=1; k<=24; k++)
    {

        int n = pow(2,k);

        double time_start = omp_get_wtime();

        double sumR = 0.0;

        #pragma omp parallel for num_threads(nthreads) reduction(+:sumR)
            for(int j=1; j<=n; j++){
                sumR += 1.0/((double)j*(double)j);
            }

        double piR = sqrt(6*sumR);
        double error = fabs(piR-pi);
        double duration = omp_get_wtime() - time_start;

        out_data1 << k << " " << scientific << setprecision(18) << log(error) << endl;
        out_data2 << k << " " << scientific << setprecision(18) << duration << endl;

    }

    return 0;



}




