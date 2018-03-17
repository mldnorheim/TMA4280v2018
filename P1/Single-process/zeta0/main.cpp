#include <iostream>
#include <math.h>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <fstream>
/*
#include "zeta0.cpp"
#include "mach0.cpp"
*/
#include "verification.h"

using namespace std;

//template <typename T>
//std::vector<T> linspace(T a, T b, size_t N) {
//    T h = (b - a) / static_cast<T>(N-1);
//    std::vector<T> xs(N);
//    typename std::vector<T>::iterator x;
//    T val;
//    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h)
//        *x = val;
//    return xs;
//}

int main(int argc, char **argv)
{
    //hello();
    if (argc < 3) { printf("Missing argument \n"); return 1; }
    printf("my program name is %s\n", argv[0]);
    int type = atoi(argv[1]);
    int n    = atoi(argv[2]);

    //cout << "Please enter a value for n: " << endl;
    //int n;
    //cin >> n;



    zeta0(n, type);

    mach0(n, type);

    return 0;
}
