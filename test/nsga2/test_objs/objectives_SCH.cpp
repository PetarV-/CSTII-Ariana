#include <functional>
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

vector<function<double(vector<double>)> > get_objectives()
{
    vector<function<double(vector<double>)> > ret;
    ret.push_back(SCH_F1);
    ret.push_back(SCH_F2);
    return ret;
}
