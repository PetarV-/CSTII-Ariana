#include <vector>
#include <math.h>

#include <objectives.h>

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

vector<pfunc_t> get_objectives()
{
    vector<pfunc_t> ret;
    ret.push_back(&FON_F1);
    ret.push_back(&FON_F2);
    return ret;
}