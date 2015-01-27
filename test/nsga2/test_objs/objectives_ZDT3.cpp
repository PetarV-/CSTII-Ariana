#include <math.h>
#include <functional>
#include <vector>

#include <objectives.h>

using namespace std;

double ZDT3_F1(vector<double> X)
{
    return X[0];
}

double ZDT3_F2(vector<double> X)
{
    double S = 0.0;
    for (int i=1;i<30;i++) S += X[i];
    double gX = 1.0 + 9.0 * S / 29.0;
    return gX * (1.0 - sqrt(X[0] / gX) - X[0] * sin(10.0 * M_PI * X[0]) / gX);
}

vector<function<double(vector<double>)> > get_objectives()
{
    vector<function<double(vector<double>)> > ret;
    ret.push_back(ZDT3_F1);
    ret.push_back(ZDT3_F2);
    return ret;
}
