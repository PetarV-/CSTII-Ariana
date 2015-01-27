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

#include <hmm.h>

#define EPS 1e-3
#define LIM 1e-20
#define INF 1e20

#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;
typedef unsigned int uint;
typedef long long lld;
typedef unsigned long long llu;

tdHMM::tdHMM(int n, int k, int obs, double *P, double **A, double **B) : HMM(n*k, obs), k(k)
{
    this -> P_init = new double[this -> k];
    for (int i=0;i<k;i++) P_init[i] = P[i];
    
    this -> dim_trans = new double*[k];
    for (int i=0;i<k;i++)
    {
        this -> dim_trans[i] = new double[k];
        for (int j=0;j<k;j++)
        {
            this -> dim_trans[i][j] = A[i][j];
        }
    }
    
    this -> dim_obs = new double*[k];
    for (int i=0;i<k;i++)
    {
        this -> dim_obs[i] = new double[obs];
        for (int j=0;j<obs;j++)
        {
            this -> dim_obs[i][j] = B[i][j];
        }
    }
    
    for (int i=0;i<n*k;i++)
    {
        for (int j=0;j<n*k;j++)
        {
            this -> T[i][j] = A[i/n][j/n];
        }
    }
    
    for (int i=0;i<n*k;i++)
    {
        for (int j=0;j<obs;j++)
        {
            this -> P[i][j] = B[i/n][j];
        }
    }
}
