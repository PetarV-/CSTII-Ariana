#include <vector>
#include <math.h>

#include "objectives.h"

using namespace std;

double ZDT2_F1(vector<double> X)
{
    return X[0];
}

double ZDT2_F2(vector<double> X)
{
    double S = 0.0;
    for (int i=1;i<30;i++) S += X[i];
    double gX = 1.0 + 9.0 * S / 29.0;
    return gX * (1.0 - (X[0] / gX) * (X[0] / gX));
}

vector<pfunc_t> get_objectives()
{
    vector<pfunc_t> ret;
    ret.push_back(&ZDT2_F1);
    ret.push_back(&ZDT2_F2);
    return ret;
}