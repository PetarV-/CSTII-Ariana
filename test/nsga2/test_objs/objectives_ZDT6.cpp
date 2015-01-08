#include <vector>
#include <math.h>

#include "objectives.h"

using namespace std;

double ZDT6_F1(vector<double> X)
{
    return 1.0 - exp(-4.0 * X[0]) * pow(sin(4.0 * M_PI * X[0]), 6);
}

double ZDT6_F2(vector<double> X)
{
    double v = ZDT6_F1(X);
    double S = 0.0;
    for (int i=1;i<10;i++) S += X[i];
    double gX = 1.0 + 9.0 * pow(S / 9.0, 0.25);
    return gX * (1.0 - (v/gX) * (v/gX));
}

vector<pfunc_t> get_objectives()
{
    vector<pfunc_t> ret;
    ret.push_back(&ZDT6_F1);
    ret.push_back(&ZDT6_F2);
    return ret;
}