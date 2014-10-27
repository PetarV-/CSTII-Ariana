#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <iostream>

#define EPS 1e-3

using namespace std;

class HMM
{
public:
    int n; // number of nodes
    int obs; // number of observations
    double **T; // transition probability matrix
    double **P; // output probability matrix
    
    HMM(int n, int obs, double **T, double **P)
    {
        this -> n = n;
        this -> obs = obs;
        
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