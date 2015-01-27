#include <math.h>
#include <functional>
#include <vector>

#include <objectives.h>

using namespace std;

double ZDT4_F1(vector<double> X)
{
    return X[0];
}

double ZDT4_F2(vector<double> X)
{
    double S = 0.0;
    for (int i=1;i<10;i++) S += X[i] * X[i] - 10.0 * cos(4.0 * M_PI * X[i]);
    double gX = 1.0 + 10.0 * 9.0 + S;
    return gX * (1.0 - sqrt(X[0] / gX));
}

vector<function<double(vector<double>)> > get_objectives()
{
    vector<function<double(vector<double>)> > ret;
    ret.push_back(ZDT4_F1);
    ret.push_back(ZDT4_F2);
    return ret;
}
