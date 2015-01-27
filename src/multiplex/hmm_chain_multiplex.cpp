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
typedef unsigned int uint;
typedef long long lld;
typedef unsigned long long llu;

HMMChainMultiplex::HMMChainMultiplex(int obs, int L) : obs(obs), L(L)
{
    this -> layers.resize(L);
    for (int i=0;i<L;i++)
    {
        this -> layers[i] = new SimpleChainGMHMM(obs);
    }
    
    this -> omega = new double*[L];
    for (int i=0;i<L;i++)
    {
        this -> omega[i] = new double[L];
        for (int j=0;j<L;j++)
        {
            this -> omega[i][j] = 0.0;
        }
    }
}

void HMMChainMultiplex::set_omega(double **omega)
{
    for (int i=0;i<L;i++)
    {
        for (int j=0;j<L;j++)
        {
            this -> omega[i][j] = omega[i][j];
        }
    }
}

void HMMChainMultiplex::train(vector<vector<vector<double> > > train_set)
{
    for (int l=0;l<L;l++)
    {
        vector<vector<double> > curr_set(train_set.size(), vector<double>(obs));
        for (uint i=0;i<train_set.size();i++)
        {
            for (int j=0;j<obs;j++)
            {
                curr_set[i][j] = train_set[i][j][l];
            }
        }
        layers[l] -> train(curr_set);
    }
    
    
}

double HMMChainMultiplex::log_likelihood(vector<vector<double> > test_data)
{
    
}
