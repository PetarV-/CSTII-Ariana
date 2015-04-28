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

HMM::HMM(int n, int obs) : n(n), obs(obs)
{
    this -> pi = new double[n];
    for (int i=0;i<n;i++)
    {
        this -> pi[i] = 1.0 / n;
    }
    
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
    
HMM::HMM(int n, int obs, double *pi, double **T, double **P) : n(n), obs(obs)
{
    double sum = 0.0;
    this -> pi = new double[n];
    for (int i=0;i<n;i++)
    {
        this -> pi[i] = pi[i];
        sum += pi[i];
    }
    assert(fabs(sum - 1.0) < EPS);
    
    this -> T = new double*[n];
    for (int i=0;i<n;i++)
    {
        this -> T[i] = new double[n];
        sum = 0.0;
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
        sum = 0.0;
        for (int j=0;j<obs;j++)
        {
            this -> P[i][j] = P[i][j];
            sum += P[i][j];
        }
        assert(fabs(sum - 1.0) < EPS);
    }
}

HMM::~HMM()
{
    for (int i=0;i<n;i++) delete[] T[i];
    delete[] T;
    
    for (int i=0;i<n;i++) delete[] P[i];
    delete[] P;
    
    delete[] pi;
}

void HMM::forward_backward(vector<int> &Y)
{
    int tlen = Y.size();
        
    // forward
    A = new double*[tlen];
    Arn = new int[tlen];
    for (int i=0;i<tlen;i++) A[i] = new double[n];
    for (int i=0;i<n;i++) A[0][i] = pi[i] * P[i][Y[0]];
    Arn[0] = 0;
    for (int t=1;t<tlen;t++)
    {
        double fullsum = 0.0;
        for (int j=0;j<n;j++)
        {
            double sum = 0.0;
            for (int i=0;i<n;i++)
            {
                sum += A[t-1][i] * T[i][j] * P[j][Y[t]];
            }
            A[t][j] = sum;
            fullsum += sum;
        }
        Arn[t] = Arn[t-1];
        if (fullsum < LIM)
        {
            // renormalise as necessary
            Arn[t]++;
            for (int j=0;j<n;j++) A[t][j] *= INF;
        }
    }
        
    // backward
    B = new double*[tlen];
    Brn = new int[tlen];
    for (int i=0;i<tlen;i++) B[i] = new double[n];
    for (int i=0;i<n;i++) B[tlen-1][i] = 1.0;
    Brn[tlen-1] = 0;
    for (int t=tlen-2;t>=0;t--)
    {
        double fullsum = 0.0;
        for (int i=0;i<n;i++)
        {
            double sum = 0.0;
            for (int j=0;j<n;j++)
            {
                sum += T[i][j] * P[j][Y[t+1]] * B[t+1][j];
            }
            B[t][i] = sum;
            fullsum += sum;
        }
        Brn[t] = Brn[t+1];
        if (fullsum < LIM)
        {
            // renormalise as necessary
            Brn[t]++;
            for (int j=0;j<n;j++) B[t][j] *= INF;
        }
    }
    
    // likelihood
    likelihood = 0.0;
    for (int i=0;i<n;i++) likelihood += A[0][i] * B[0][i];
    renorm = Arn[0] + Brn[0];
    while (likelihood < LIM)
    {
        likelihood *= INF;
        renorm++;
    }
    
    Pst = new double*[tlen];
    for (int i=0;i<tlen;i++) Pst[i] = new double[n];
    for (int t=0;t<tlen;t++)
    {
        double sum = 0.0;
        for (int i=0;i<n;i++)
        {
            Pst[t][i] = A[t][i] * B[t][i];
            sum += Pst[t][i];
        }
        for (int i=0;i<n;i++) Pst[t][i] /= sum;
    }
}

double** HMM::get_A()
{
    return A;
}

double** HMM::get_B()
{
    return B;
}

double* HMM::get_pi()
{
    return pi;
}

double** HMM::get_T()
{
    return T;
}

double** HMM::get_P()
{
    return P;
}

vector<int> HMM::viterbi(vector<int> &Y)
{
    int tlen = Y.size();
        
    double **V = new double*[tlen];
    for (int i=0;i<tlen;i++) V[i] = new double[n];
    
    int **prev = new int*[tlen];
    for (int i=0;i<tlen;i++) prev[i] = new int[n];
    
    // initialisation
    for (int i=0;i<n;i++)
    {
        V[0][i] = log(pi[i]) + log(P[i][Y[0]]);
        prev[0][i] = 0;
    }
    
    for (int t=1;t<tlen;t++)
    {
        for (int i=0;i<n;i++)
        {
            double maxx = -1.0;
            int maxState = -1;
            for (int j=0;j<n;j++)
            {
                double curr = V[t-1][j] + log(T[j][i]) + log(P[i][Y[t]]);
                if (maxState == -1 || curr > maxx)
                {
                    maxx = curr;
                    maxState = j;
                }
            }
            V[t][i] = maxx;
            prev[t][i] = maxState;
        }
    }
        
    double best = -1;
    int bestState = -1;
    for (int i=0;i<n;i++)
    {
        if (bestState == -1 || V[tlen-1][i] > best)
        {
            best = V[tlen-1][i];
            bestState = i;
        }
    }
    
    vector<int> ret;
    ret.resize(tlen);
    for (int t=tlen-1;t>=0;t--)
    {
        ret[t] = bestState;
        bestState = prev[t][bestState];
    }
    
    for (int i=0;i<tlen;i++) delete[] V[i];
    delete[] V;
    
    for (int i=0;i<tlen;i++) delete[] prev[i];
    delete[] prev;
    
    return ret;
}

void HMM::baumwelch(vector<int> &Y)
{
    forward_backward(Y);
        
    int tlen = Y.size();
        
    double **nextP = new double*[n];
    for (int i=0;i<n;i++) nextP[i] = new double[obs];
        
    double powers[10];
    for (int i=0;i<10;i++) powers[i] = pow(LIM, i-6);
        
    double PP, QQ;
    
    for (int i=0;i<n;i++)
    {
        pi[i] = ((A[0][i] * B[0][i]) / likelihood) * powers[Arn[0] + Brn[0] - renorm + 6];
        QQ = 0.0;
        for (int k=0;k<obs;k++)
        {
            nextP[i][k] = 0.0;
        }
        for (int t=0;t<tlen-1;t++)
        {
            double curr = ((A[t][i] * B[t][i]) / likelihood) * powers[Arn[t] + Brn[t] - renorm + 6];
            QQ += curr;
            nextP[i][Y[t]] += curr;
        }
        for (int j=0;j<n;j++)
        {
            PP = 0.0;
            for (int t=0;t<tlen-1;t++)
            {
                PP += A[t][i] * P[j][Y[t+1]] * B[t+1][j] * powers[Arn[t] + Brn[t] - renorm + 6] / likelihood;
            }
            T[i][j] *= PP / QQ;
        }
        for (int k=0;k<obs;k++)
        {
            nextP[i][k] /= QQ;
        }
    }
    
    for (int i=0;i<n;i++) delete[] P[i];
    delete[] P;
    
    P = nextP;
}
