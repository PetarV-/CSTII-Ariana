/*
 Petar 'PetarV' Velickovic
 Data Structure: Hidden Markov Model
*/

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

#define EPS 1e-3

#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;
typedef long long lld;
typedef unsigned long long llu;

class HMM
{
public:
    int n; // number of nodes
    int obs; // number of observations
    double **T; // transition probability matrix
    double **P; // output probability matrix
    
    double **A, **B; // alpha, beta matrices
    
    
    HMM(int n, int obs) : n(n), obs(obs)
    {
        this -> T = new double*[n];
        for (int i=0;i<n;i++)
        {
            this -> T[i] = new double[n];
            for (int j=0;j<n;j++)
            {
                this -> T[i][j] = 1.0 / n;
            }
        }
        
        this -> P = new double*[n];
        for (int i=0;i<n;i++)
        {
            this -> P[i] = new double[obs];
            for (int j=0;j<obs;j++)
            {
                this -> T[i][j] = 1.0 / obs;
            }
        }
    }
    
    HMM(int n, int obs, double **T, double **P) : n(n), obs(obs)
    {
        this -> T = new double*[n];
        for (int i=0;i<n;i++)
        {
            this -> T[i] = new double[n];
            double sum = 0.0;
            for (int j=0;j<n;j++)
            {
                this -> T[i][j] = T[i][j];
                sum += T[i][j];
            }
            assert(fabs(sum - 1.0) < EPS);
        }
        
        this -> P = new double*[n];
        for (int i=0;i<n;i++)
        {
            this -> P[i] = new double[obs];
            double sum = 0.0;
            for (int j=0;j<obs;j++)
            {
                this -> P[i][j] = P[i][j];
                sum += P[i][j];
            }
            assert(fabs(sum - 1.0) < EPS);
        }
    }
    
    void forward_backward(vector<int> Y)
    {
        int tlen = Y.size();
        
        A = new double*[tlen];
        for (int i=0;i<tlen;i++) A[i] = new double[n];
        for (int i=0;i<n;i++) A[0][i] = P[i][Y[0]];
        for (int t=1;t<tlen;t++)
        {
            for (int j=0;j<n;j++)
            {
                double sum = 0.0;
                for (int i=0;i<n;i++)
                {
                    sum += A[t-1][i] * T[i][j] * P[j][Y[t]];
                }
                A[t][j] = sum;
            }
        }
        
        B = new double*[tlen];
        for (int i=0;i<tlen;i++) B[i] = new double[n];
        for (int i=0;i<n;i++) B[tlen-1][i] = 1.0;
        for (int t=tlen-2;t>=0;t--)
        {
            for (int i=0;i<n;i++)
            {
                double sum = 0.0;
                for (int j=0;j<n;j++)
                {
                    sum += T[i][j] * P[j][Y[t+1]] * B[t+1][j];
                }
                B[t][i] = sum;
            }
        }
    }
};

int main()
{
    // test
    int n = 2;
    int obs = 2;
    double **T = new double*[n];
    for (int i=0;i<n;i++)
    {
        T[i] = new double[n];
        for (int j=0;j<n;j++)
        {
            T[i][j] = 0.5;
        }
    }
    double **P = new double*[n];
    for (int i=0;i<n;i++)
    {
        P[i] = new double[obs];
        for (int j=0;j<obs;j++)
        {
            P[i][j] = 0.5;
        }
    }
    HMM *x = new HMM(n, obs, T, P);
}