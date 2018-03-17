#include "verification.h"
#include <fstream>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>
#include <cmath>


using namespace std;


void mach0(int n, int type)      // MACHIN:
{

  //  int type;
  //  int n;


    double SnM5 = 0.0;
    double SnM239 = 0.0;
    for(int j=1; j<=n; j++){
        // x=1/5
        double vjx5 = pow(-1.0,j-1) * pow(1.0/5.0,2.0*j-1) / (2.0*j-1.0);
        //cout << vjx5 << endl;
        SnM5 += vjx5;

        // x=1/239
        double vjx239 = pow(-1.0,j-1) * pow(1.0/239.0,2.0*j-1) / (2.0*j-1.0);
        //cout << vjx239 << endl << "---" << endl;
        SnM239 += vjx239;
    }
    double piM = 4.0*(4.0*SnM5 - SnM239);


    // TEST
    if (type == 0)
    {
        cout << "Machin approximates pi as: " << fixed << scientific << setprecision(15) << piM << endl;
    }


    // UNIT TEST, n=3:
    else if (type == 1)
    {
        double utpiM = 38279241713339684.0/12184551018734375.0;
        double diffM = piM - utpiM;
        cout << "Machin unit test with n=3 gives approximation: " << scientific << setprecision(15) << utpiM << endl;
        cout << "and the difference is: " << fixed << scientific << setprecision(15) << diffM << endl;
    }


    // VERIFICATION TEST, n=2^k for k=1,...,24:
    else if (type == 2)
    {
        ofstream out_data("vtestmachlog.txt");

        const double pi = 3.141592653589793;         // "exact" pi
        //vector<double> PIM (24,0);
        double PIM = 0.0;
        vector<double> errorM (24,0.0);
        for(int k=1; k<=24; k++){
            int N = pow(2,k);
            double SNM5 = 0.0;
            double SNM239 = 0.0;
            for(int j=1; j<=N; j++){
                // x=1/5
                double Vjx5 = pow(-1.0,j-1) * pow(1.0/5.0,2.0*j-1) / (2.0*j-1.0);
                SNM5 += Vjx5;
                // x=1/239
                double Vjx239 = pow(-1.0,j-1) * pow(1.0/239.0,2.0*j-1) / (2.0*j-1.0);
                SNM239 += Vjx239;
            }
            PIM = 4.0*(4.0*SNM5 - SNM239);
            errorM[k] = fabs(pi-PIM);
            out_data << k << " " << scientific << setprecision(18) << log(fabs(pi-PIM)) << endl;
            //cout << "n = " << N << ", |pi-pi_n_mach| = " << abs(pi-PIM) << endl;
        }
        //cout << "The errors for Machin and different N's are: " << errorM << endl;
    }

}
