#include <math.h>
#include <functional>
#include <vector>

using namespace std;

double FON_F1(vector<double> X)
{
    double sum = 0.0;
    for (int i=0;i<3;i++)
    {
        sum += (X[i] - 1.0/sqrt(3.0)) * (X[i] - 1.0 / sqrt(3.0));
    }
    return 1.0 - exp(-sum);
}

double FON_F2(vector<double> X)
{
    double sum = 0.0;
    for (int i=0;i<3;i++)
    {
        sum += (X[i] + 1.0/sqrt(3.0)) * (X[i] + 1.0 / sqrt(3.0));
    }
    return 1.0 - exp(-sum);
}

vector<function<double(vector<double>)> > get_objectives()
{
    vector<function<double(vector<double>)> > ret;
    ret.push_back(FON_F1);
    ret.push_back(FON_F2);
    return ret;
}
