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
#define LIM 1e-20
#define INF 1e20

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
    
    double **A, **B, **Pst; // alpha, beta, state prob. matrices
    int *Arn, *Brn; // how many times have we renormalised A, B
    int renorm;
    double likelihood;
    
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
        
        // forward
        A = new double*[tlen];
        Arn = new int[tlen];
        for (int i=0;i<tlen;i++) A[i] = new double[n];
        for (int i=0;i<n;i++) A[0][i] = P[i][Y[0]];
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
    
    vector<int> viterbi(vector<int> Y)
    {
        int tlen = Y.size();
        
        double **V = new double*[tlen];
        for (int i=0;i<tlen;i++) V[i] = new double[n];
        int ***paths = new int**[2];
        for (int id=0;id<2;id++)
        {
            paths[id] = new int*[n];
            for (int i=0;i<n;i++) paths[id][i] = new int[tlen];
        }
        
        // initialisation
        for (int i=0;i<n;i++)
        {
            V[0][i] = P[i][Y[0]];
            paths[0][i][0] = i;
        }
        
        for (int t=1;t<tlen;t++)
        {
            for (int i=0;i<n;i++)
            {
                double maxx = -1.0;
                int maxState = -1;
                for (int j=0;j<n;j++)
                {
                    double curr = V[t-1][j] * T[j][i] * P[i][Y[t]];
                    if (curr > maxx)
                    {
                        maxx = curr;
                        maxState = j;
                    }
                }
                V[t][i] = maxx;
                for (int j=0;j<t;j++)
                {
                    paths[t&1][i][j] = paths[(t+1)&1][maxState][j];
                }
                paths[t&1][i][t] = i;
            }
        }
        
        double best = -1;
        int bestState = -1;
        for (int i=0;i<n;i++)
        {
            if (V[tlen-1][i] > best)
            {
                best = V[tlen-1][i];
                bestState = i;
            }
        }
        
        vector<int> ret;
        ret.resize(tlen);
        for (int i=0;i<tlen;i++) ret[i] = paths[(tlen-1)&1][bestState][i];
        
        return ret;
    }
    
    void baumwelch(vector<int> Y)
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
};

int main()
{
    // test
    int n = 2;
    int obs = 3;
    
    double **T = new double*[n];
    for (int i=0;i<n;i++)
    {
        T[i] = new double[n];
    }
    T[0][0] = 0.7; T[0][1] = 0.3;
    T[1][0] = 0.4; T[1][1] = 0.6;
    
    double **P = new double*[n];
    for (int i=0;i<n;i++)
    {
        P[i] = new double[obs];
    }
    P[0][0] = 0.5; P[0][1] = 0.4; P[0][2] = 0.1;
    P[1][0] = 0.1; P[1][1] = 0.3; P[1][2] = 0.6;
    
    HMM *x = new HMM(n, obs, T, P);
    
    vector<int> observations;
    observations.push_back(0);
    observations.push_back(1);
    observations.push_back(2);
    
    vector<int> ret = x -> viterbi(observations);
    
    for (int i=0;i<ret.size();i++) printf("%d ", ret[i]);
    printf("\n");
}