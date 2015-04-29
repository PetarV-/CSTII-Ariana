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
#include <tuple>

#include <hmm.h>

#define EPS 1e-3

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
   
    printf("Testing forward/backward algorithm... ");
    
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
    
    double exp_A[5][2];
    exp_A[0][0] = 0.8182; exp_A[0][1] = 0.1818;
    exp_A[1][0] = 0.8834; exp_A[1][1] = 0.1166;
    exp_A[2][0] = 0.1907; exp_A[2][1] = 0.8093;
    exp_A[3][0] = 0.7308; exp_A[3][1] = 0.2692;
    exp_A[4][0] = 0.8673; exp_A[4][1] = 0.1327;
    
    double exp_B[5][2];
    exp_B[0][0] = 0.5923; exp_B[0][1] = 0.4077;
    exp_B[1][0] = 0.3763; exp_B[1][1] = 0.6237;
    exp_B[2][0] = 0.6533; exp_B[2][1] = 0.3467;
    exp_B[3][0] = 0.6273; exp_B[3][1] = 0.3727;
    exp_B[4][0] = 0.5000; exp_B[4][1] = 0.5000;
    
    HMM *x = new HMM(n, obs, pi, T, P);
    
    vector<int> observations;
    observations.push_back(0);
    observations.push_back(0);
    observations.push_back(1);
    observations.push_back(0);
    observations.push_back(0);
    
    //printf("Running forward-backward algorithm... ");
    tuple<double**, double*, double> ret = x -> forward(observations);
    double **AA = get<0>(ret);
    double **BB = x -> backward(observations, get<1>(ret));
    
    //printf("Forward likelihoods:\n");
    for (int t=0;t<5;t++)
    {
        double sum = 0.0;
        for (int i=0;i<n;i++)
        {
            sum += AA[t][i];
        }
        for (int i=0;i<n;i++)
        {
            assert(fabs(exp_A[t][i] - AA[t][i] / sum) < EPS);
            //printf("%lf ", AA[t][i] / sum);
        }
        //printf("\n");
    }
    
    //printf("Backward likelihoods:\n");
    for (int t=0;t<5;t++)
    {
        double sum = 0.0;
        for (int i=0;i<n;i++)
        {
            sum += BB[t][i];
        }
        for (int i=0;i<n;i++)
        {
            assert(fabs(exp_B[t][i] - BB[t][i] / sum) < EPS);
            //printf("%lf ", BB[t][i] / sum);
        }
        //printf("\n");
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
    
    printf("\033[1;32mOK!\033[0m\n");
    
    printf("Testing Viterbi algorithm... ");
    
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
    
    int exp_seq[3];
    exp_seq[0] = 0; exp_seq[1] = 0; exp_seq[2] = 1;
    
    x = new HMM(n, obs, pi, T, P);
    
    observations.push_back(0);
    observations.push_back(1);
    observations.push_back(2);
    
    //printf("Running Viterbi algorithm...\n");
    vector<int> seq = x -> viterbi(observations);
    
    //printf("Most likely state sequence: ");
    for (uint i=0;i<seq.size();i++)
    {
        assert(seq[i] == exp_seq[i]);
        //printf("%d ", ret[i]);
    }
    //printf("\n");
    
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
    
    printf("\033[1;32mOK!\033[0m\n");
    
    printf("Testing Baum-Welch algorithm... ");
    
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
    
    double exp_pi[2];
    exp_pi[0] = 1.0000; exp_pi[1] = 0.0000;
    
    double exp_T[2][2];
    exp_T[0][0] = 0.6906; exp_T[0][1] = 0.3091;
    exp_T[1][0] = 0.0934; exp_T[1][1] = 0.9066;
    
    double exp_P[2][3];
    exp_P[0][0] = 0.5807; exp_P[0][1] = 0.0010; exp_P[0][2] = 0.4183;
    exp_P[1][0] = 0.0000; exp_P[1][1] = 0.7621; exp_P[1][2] = 0.2379;
    
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
    
    vector<vector<int> > obs1;
    obs1.push_back(observations);
    
    //printf("Running Baum-Welch algorithm...\n");
    x -> baumwelch(obs1, 500, 1e-10);
    
    //printf("Start-state probability vector:\n");
    pi = x -> get_pi();
    for (int i=0;i<n;i++)
    {
        //printf("%lf ", pi[i]);
        assert(fabs(pi[i] - exp_pi[i]) < EPS);
    }
    //printf("\n");
    
    //printf("Transition probability matrix:\n");
    T = x -> get_T();
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            //printf("%lf ", T[i][j]);
            assert(fabs(T[i][j] - exp_T[i][j]) < EPS);
        }
        //printf("\n");
    }
    
    //printf("Emission probability matrix:\n");
    P = x -> get_P();
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<obs;j++)
        {
            //printf("%lf ", P[i][j]);
            assert(fabs(P[i][j] - exp_P[i][j]) < EPS);
        }
        //printf("\n");
    }
    
    printf("\033[1;32mOK!\033[0m\n");
    
    printf("\033[1;32mSUCCESS!!!\033[0m\n");
    
    return 0;
}
