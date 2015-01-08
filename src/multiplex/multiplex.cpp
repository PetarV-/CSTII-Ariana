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

#include "multiplex.h"

#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;
typedef long long lld;
typedef unsigned long long llu;

Multiplex::Multiplex(int n, vector<AbstractGraphLayer> lyrs)
{
    this -> n = n;
    for (int i=0;i<lyrs.size();i++) assert(lyrs[i].get_n() == n);
    M = new double***[lyrs.size()];
    for (int i=0;i<lyrs.size();i++)
    {
        M[i] = new double**[lyrs.size()];
        for (int j=0;j<lyrs.size();j++)
        {
            M[i][j] = new double*[n];
            for (int k=0;k<n;k++)
            {
                M[i][j][k] = new double[n];
                for (int l=0;l<n;l++)
                {
                    if (i != j) M[i][j][k][l] = 0.0;
                    else M[i][j][k][l] = lyrs[i].get_adj(k, l);
                }
            }
        }
    }
}

double Multiplex::get_edge(int lyr1, int lyr2, int n1, int n2)
{
    return M[lyr1][lyr2][n1][n2];
}

void Multiplex::set_edge(int lyr1, int lyr2, int n1, int n2, double val)
{
    M[lyr1][lyr2][n1][n2] = val;
}

int main()
{
    
}