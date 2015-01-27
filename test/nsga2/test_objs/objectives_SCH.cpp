#include <vector>

#include <objectives.h>

using namespace std;

double SCH_F1(vector<double> X)
{
    return X[0] * X[0];
}

double SCH_F2(vector<double> X)
{
    return (X[0] - 2) * (X[0] - 2);
}

vector<pfunc_t> get_objectives()
{
    vector<pfunc_t> ret;
    ret.push_back(&SCH_F1);
    ret.push_back(&SCH_F2);
    return ret;
}