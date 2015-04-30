#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <deque>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <complex>

#include <layers.h>

#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;
typedef unsigned int uint;
typedef long long lld;
typedef unsigned long long llu;

SimpleGraphLayer::SimpleGraphLayer(int n) : n(n)
{
    this -> G = new double*[n];
    for (int i=0;i<n;i++)
    {
        G[i] = new double[n];
        for (int j=0;j<n;j++)
        {
            if (i == j) G[i][j] = 1.0;
            else G[i][j] = 0.0;
        }
    }
}

SimpleGraphLayer::SimpleGraphLayer(int n, double **A) : n(n)
{
    this -> G = new double*[n];
    for (int i=0;i<n;i++)
    {
        G[i] = new double[n];
        double sum = 0.0;
        for (int j=0;j<n;j++)
        {
            G[i][j] = A[i][j];
            sum += A[i][j];
        }
        for (int j=0;j<n;j++)
        {
            G[i][j] /= A[i][j];
        }
    }
}

SimpleGraphLayer::~SimpleGraphLayer()
{
    for (int i=0;i<n;i++) delete[] G[i];
    delete[] G;
}

int SimpleGraphLayer::get_n()
{
    return n;
}

double SimpleGraphLayer::get_adj(int x, int y)
{
    return G[x][y];
}

void SimpleGraphLayer::set_adj(int x, int y, double val)
{
    G[x][y] = val;
}

void SimpleGraphLayer::train(vector<vector<double> > &train_set)
{
    // Take the mean of all distances, then normalise that
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            G[i][j] = 0.0;
        }
    }
    
    double epsilon = 1e-6; // the smallest score possible, can be tweaked
    double max_diff = 0.0; // the maximal difference observed
    
    for (uint i=0;i<train_set.size();i++)
    {
        for (int j=0;j<n;j++)
        {
            for (int k=0;k<n;k++)
            {
                double curr_diff = fabs(train_set[i][j] - train_set[i][k]);
                if (curr_diff > max_diff) max_diff = curr_diff;
            }
        }
    }
    
    double lambda = -log(epsilon) / max_diff;
    
    for (uint i=0;i<train_set.size();i++)
    {
        for (int j=0;j<n;j++)
        {
            for (int k=0;k<n;k++)
            {
                if (j != k) G[j][k] += exp(-lambda * fabs(train_set[i][j] - train_set[i][k])); // e^(-lambda*diff)
            }
        }
    }
    
    for (int i=0;i<n;i++)
    {
        double sum = 0.0;
        for (int j=0;j<n;j++)
        {
            sum += G[i][j];
        }
        for (int j=0;j<n;j++)
        {
            G[i][j] /= sum;
        }
    }
}

double SimpleGraphLayer::log_likelihood(vector<double> &test_data)
{
    // first build another simple layer using the test data
    vector<vector<double> > test_set(1, test_data);
    SimpleGraphLayer *test_lyr = new SimpleGraphLayer(n);
    
    test_lyr -> train(test_set);
    
    // now compare
    double ret = 0.0;
    
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            ret += (this -> G[i][j] - test_lyr -> G[i][j]) * (this -> G[i][j] - test_lyr -> G[i][j]);
        }
    }
    
    delete test_lyr;
    
    // as the similarity is required rather than distance,
    // return the negative log
    return -log(sqrt(ret));
}
