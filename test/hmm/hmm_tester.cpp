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

#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;
typedef long long lld;
typedef unsigned long long llu;

int main()
{
    printf("Running HMM Tests...\n");
   
    printf("Testing forward/backward algorithm...\n");
    
    int n = 2;
    int obs = 2;
    
    double *pi = new double[n];
    for (int i=0;i<n;i++)
    {
        pi[i] = 0.5;
    }
    
    double **T = new double*[n];
    for (int i=0;i<n;i++)
    {
        T[i] = new double[n];
    }
    T[0][0] = 0.7; T[0][1] = 0.3;
    T[1][0] = 0.3; T[1][1] = 0.7;
    
    double **P = new double*[n];
    for (int i=0;i<n;i++)
    {
        P[i] = new double[obs];
    }
    P[0][0] = 0.9; P[0][1] = 0.1;
    P[1][0] = 0.2; P[1][1] = 0.8;
    
    HMM *x = new HMM(n, obs, pi, T, P);
    
    vector<int> observations;
    observations.push_back(0);
    observations.push_back(0);
    observations.push_back(1);
    observations.push_back(0);
    observations.push_back(0);
    
    printf("Running forward-backward algorithm...\n");
    x -> forward_backward(observations);
    
    printf("Forward likelihoods:\n");
    double **AA = x -> get_A();
    double **BB = x -> get_B();
    for (int t=0;t<5;t++)
    {
        double sum = 0.0;
        for (int i=0;i<n;i++)
        {
            sum += AA[t][i];
        }
        for (int i=0;i<n;i++)
        {
            printf("%lf ", AA[t][i] / sum);
        }
        printf("\n");
    }
    
    printf("Backward likelihoods:\n");
    for (int t=0;t<5;t++)
    {
        double sum = 0.0;
        for (int i=0;i<n;i++)
        {
            sum += BB[t][i];
        }
        for (int i=0;i<n;i++)
        {
            printf("%lf ", BB[t][i] / sum);
        }
        printf("\n");
    }
    
    delete x;
    for (int i=0;i<n;i++)
    {
        delete[] T[i];
        delete[] P[i];
    }
    delete[] pi;
    delete[] T;
    delete[] P;
    observations.clear();
    
    printf("Forward-backward testing completed!\n");
    
    printf("Testing Viterbi algorithm...\n");
    
    n = 2;
    obs = 3;
    
    pi = new double[n];
    pi[0] = 0.6; pi[1] = 0.4;
    
    T = new double*[n];
    for (int i=0;i<n;i++)
    {
        T[i] = new double[n];
    }
    T[0][0] = 0.7; T[0][1] = 0.3;
    T[1][0] = 0.4; T[1][1] = 0.6;
    
    P = new double*[n];
    for (int i=0;i<n;i++)
    {
        P[i] = new double[obs];
    }
    P[0][0] = 0.5; P[0][1] = 0.4; P[0][2] = 0.1;
    P[1][0] = 0.1; P[1][1] = 0.3; P[1][2] = 0.6;
    
    x = new HMM(n, obs, pi, T, P);
    
    observations.push_back(0);
    observations.push_back(1);
    observations.push_back(2);
    
    printf("Running Viterbi algorithm...\n");
    vector<int> ret = x -> viterbi(observations);
    
    printf("Most likely state sequence: ");
    for (uint i=0;i<ret.size();i++) printf("%d ", ret[i]);
    printf("\n");
    
    delete x;
    for (int i=0;i<n;i++)
    {
        delete[] T[i];
        delete[] P[i];
    }
    delete[] pi;
    delete[] T;
    delete[] P;
    observations.clear();
    
    printf("Viterbi testing completed!\n");
    
    printf("Testing Baum-Welch algorithm...\n");
    
    n = 2;
    obs = 3;
    
    pi = new double[n];
    for (int i=0;i<n;i++)
    {
        pi[i] = 0.5;
    }
    
    T = new double*[n];
    for (int i=0;i<n;i++)
    {
        T[i] = new double[n];
    }
    T[0][0] = 0.5; T[0][1] = 0.5;
    T[1][0] = 0.5; T[1][1] = 0.5;
    
    P = new double*[n];
    for (int i=0;i<n;i++)
    {
        P[i] = new double[obs];
    }
    P[0][0] = 0.4; P[0][1] = 0.1; P[0][2] = 0.5;
    P[1][0] = 0.1; P[1][1] = 0.5; P[1][2] = 0.4;
    
    x = new HMM(n, obs, pi, T, P);
    
    observations.push_back(2);
    observations.push_back(0);
    observations.push_back(0);
    observations.push_back(2);
    observations.push_back(1);
    observations.push_back(2);
    observations.push_back(1);
    observations.push_back(1);
    observations.push_back(1);
    observations.push_back(2);
    observations.push_back(1);
    observations.push_back(1);
    observations.push_back(1);
    observations.push_back(1);
    observations.push_back(1);
    observations.push_back(2);
    observations.push_back(2);
    observations.push_back(0);
    observations.push_back(0);
    observations.push_back(1);
    
    printf("Running Baum-Welch algorithm...\n");
    for (int i=0;i<100;i++) x -> baumwelch(observations);
    
    printf("Start-state probability vector:\n");
    pi = x -> get_pi();
    for (int i=0;i<n;i++)
    {
        printf("%lf ", pi[i]);
    }
    printf("\n");
    
    printf("Transition probability matrix:\n");
    T = x -> get_T();
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            printf("%lf ", T[i][j]);
        }
        printf("\n");
    }
    
    printf("Emission probability matrix:\n");
    P = x -> get_P();
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<obs;j++)
        {
            printf("%lf ", P[i][j]);
        }
        printf("\n");
    }
    
    printf("Baum-Welch testing completed!\n");
    
    
    return 0;
}
