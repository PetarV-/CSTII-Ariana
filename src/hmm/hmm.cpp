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

tuple<double**, double*, double> HMM::forward(vector<int> &Y)
{
    int Ti = Y.size();
    
    double **alpha = new double*[Ti];
    for (int i=0;i<Ti;i++)
    {
        alpha[i] = new double[n];
    }
    double *c = new double[Ti];
    
    double sum = 0.0;
    for (int i=0;i<n;i++)
    {
        alpha[0][i] = pi[i] * P[i][Y[0]];
        sum += alpha[0][i];
    }
    c[0] = 1.0 / sum;
    for (int i=0;i<n;i++)
    {
        alpha[0][i] /= sum;
    }
    
    for (int t=1;t<Ti;t++)
    {
        sum = 0.0;
        for (int i=0;i<n;i++)
        {
            alpha[t][i] = 0.0;
            for (int j=0;j<n;j++)
            {
                alpha[t][i] += alpha[t-1][j] * T[j][i];
            }
            alpha[t][i] *= P[i][Y[t]];
            sum += alpha[t][i];
        }
        
        c[t] = 1.0 / sum;
        for (int i=0;i<n;i++)
        {
            alpha[t][i] /= sum;
        }
    }
    
    double log_L = 0.0;
    for (int i=0;i<Ti;i++) log_L -= log(c[i]);
    
    return make_tuple(alpha, c, log_L);
}

double** HMM::backward(vector<int> &Y, double *c)
{
    int Ti = Y.size();
    
    double **beta = new double*[Ti];
    for (int i=0;i<Ti;i++)
    {
        beta[i] = new double[n];
    }
    for (int i=0;i<n;i++) beta[Ti-1][i] = 1.0;
    
    for (int t=Ti-2;t>=0;t--)
    {
        for (int i=0;i<n;i++)
        {
            beta[t][i] = 0.0;
            for (int j=0;j<n;j++)
            {
                beta[t][i] += T[i][j] * P[j][Y[t+1]] * beta[t+1][j];
            }
            beta[t][i] *= c[t+1];
        }
    }
    
    return beta;
}

void HMM::baumwelch(vector<vector<int> > &Ys, int iterations, double tolerance)
{
    double ***alpha = new double**[Ys.size()];
    double ***beta = new double**[Ys.size()];
    double **c = new double*[Ys.size()];
    
    double PP, QQ;
    
    double lhood = 0.0;
    double oldlhood = 0.0;
    
    for (int iter=0;iter<iterations;iter++)
    {
        lhood = 0.0;
        
        for (uint l=0;l<Ys.size();l++)
        {
            tuple<double**, double*, double> x = forward(Ys[l]);
            alpha[l] = get<0>(x);
            c[l] = get<1>(x);
            lhood += get<2>(x);
            beta[l] = backward(Ys[l], c[l]);
        }
        
        double **nextO = new double*[n];
        for (int i=0;i<n;i++) nextO[i] = new double[obs];
        
        for (int i=0;i<n;i++)
        {
            pi[i] = 0.0;
            for (uint l=0;l<Ys.size();l++)
            {
                pi[i] += alpha[l][0][i] * beta[l][0][i];
            }
            pi[i] /= Ys.size();
            
            QQ = 0.0;
            
            for (int k=0;k<obs;k++)
            {
                nextO[i][k] = 0.0;
            }
            
            for (uint l=0;l<Ys.size();l++)
            {
                for (uint t=0;t<Ys[l].size()-1;t++)
                {
                    double curr = alpha[l][t][i] * beta[l][t][i];
                    QQ += curr;
                    nextO[i][Ys[l][t]] += curr;
                }
            }
            
            for (int j=0;j<n;j++)
            {
                PP = 0.0;
                for (uint l=0;l<Ys.size();l++)
                {
                    for (uint t=0;t<Ys[l].size()-1;t++)
                    {
                        PP += alpha[l][t][i] * P[j][Ys[l][t+1]] * beta[l][t+1][j] * c[l][t+1];
                    }
                }
                T[i][j] *= PP / QQ;
            }
            
            for (uint l=0;l<Ys.size();l++)
            {
                int lim = Ys[l].size() - 1;
                double curr = alpha[l][lim][i] * beta[l][lim][i];
                QQ += curr;
                nextO[i][Ys[l][lim]] += curr;
            }
            
            for (int k=0;k<obs;k++)
            {
                nextO[i][k] /= QQ;
            }
        }
        
        for (uint l=0;l<Ys.size();l++)
        {
            for (uint t=0;t<Ys[l].size();t++)
            {
                delete[] alpha[l][t];
                delete[] beta[l][t];
            }
            delete[] alpha[l];
            delete[] beta[l];
            delete[] c[l];
        }
        
        for (int i=0;i<n;i++)
        {
            delete[] P[i];
        }
        delete[] P;
        
        P = nextO;
        
        if (fabs(lhood - oldlhood) < tolerance)
        {
            break;
        }
        
        oldlhood = lhood;
    }
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
