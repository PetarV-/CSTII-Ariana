#include <vector>
#include <math.h>

#include <objectives.h>

using namespace std;

double POL_F1(vector<double> X)
{
    double A1 = 0.5 * sin(1.0) - 2.0 * cos(1.0) + sin(2.0) - 1.5 * cos(2.0);
    double A2 = 1.5 * sin(1.0) - cos(1.0) + 2.0 * sin(2.0) - 0.5 * cos(2.0);
    double B1 = 0.5 * sin(X[0]) - 2.0 * cos(X[0]) + sin(X[1]) - 1.5 * cos(X[1]);
    double B2 = 1.5 * sin(X[0]) - cos(X[0]) + 2.0 * sin(X[1]) - 0.5 * cos(X[1]);
    return 1.0 + (A1 - B1) * (A1 - B1) + (A2 - B2) * (A2 - B2);
}

double POL_F2(vector<double> X)
{
    return (X[0] + 3.0) * (X[0] + 3.0) + (X[1] + 1.0) * (X[1] + 1.0);
}

vector<pfunc_t> get_objectives()
{
    vector<pfunc_t> ret;
    ret.push_back(&POL_F1);
    ret.push_back(&POL_F2);
    return ret;
}