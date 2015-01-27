#include <vector>

using namespace std;
typedef unsigned int uint;

bool compare_euclidean(pair<vector<double>, int> &A, pair<vector<double>, int> &B)
{
    double sA = 0.0, sB = 0.0;
    for (uint i=0;i<A.first.size();i++) sA += A.first[i] * A.first[i];
    for (uint i=0;i<B.first.size();i++) sB += B.first[i] * B.first[i];
    return sA < sB;
}
