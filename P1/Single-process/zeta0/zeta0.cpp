#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fstream>
//#include "verification.h"

using namespace std;


void zeta0(int n, int type)   // RIEMANN ZETA:
{
  //  int type;
  //  int n;


    double SnR = 0.0;
    /*
    cout << "Please enter a value for n: " << endl;
    int n;
    cin >> n;
    */
    for(int j=1; j<=n; j++){
        double vj = pow(j,-2.0);
        SnR += vj;
    }
    double piR = sqrt(6.0*SnR);



    // TEST:
    if (type == 0)
    {
        cout << "Riemann approximates pi as: " << scientific << piR << endl;
    }


    // UNIT TEST, n=3:
    else if (type == 1)
    {
        double utpiR = 7.0/sqrt(6.0);
        double diffR = piR - utpiR;
        cout << "Rieman unit test with n=3 gives approximation: " << scientific << utpiR << endl;
        cout << "and the difference is: " << scientific << diffR << endl;
    }


    // VERIFICATION TEST, n=2^k for k=1,...,24:
    else if (type == 2)
    {
        ofstream out_data("vtestzetalog.txt");

        const double pi = 3.141592653589793;         // "exact" pi
        //vector<double> PIR (24,0);
        vector<double> errorR (24,0.0);
        vector<double> n (24, 0.0);
        double PIR = 0.0;
        for(int k=1; k<=24; k++){
            int N = pow(2,k);
            double SNR = 0.0;
            for(int j=1; j<=N; j++){
                double Vj = pow(j,-2.0);
                SNR += Vj;
                }
            PIR = sqrt(6.0*SNR);
            n[k] = N;
            errorR[k] = abs(pi-PIR);
            out_data << k << " " << scientific << setprecision(18) << log(fabs(pi-PIR)) << endl;
            //cout << "n = " << N << ", |pi-pi_n_zeta| = " << abs(pi-PIR) << endl;
            }
        //cout << "The errors for Riemann and different N's are: " << errorR << endl;
    }

}
