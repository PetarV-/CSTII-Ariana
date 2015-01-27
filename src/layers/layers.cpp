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
            G[i][j] = 0.0;
        }
    }
}

SimpleGraphLayer::SimpleGraphLayer(int n, double **A) : n(n)
{
    this -> G = new double*[n];
    for (int i=0;i<n;i++)
    {
        G[i] = new double[n];
        for (int j=0;j<n;j++)
        {
            G[i][j] = A[i][j];
        }
    }
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
