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
typedef unsigned int uint;
typedef long long lld;
typedef unsigned long long llu;

SimpleChainGMHMM::SimpleChainGMHMM(int obs) : obs(obs)
{
    this -> G = new double*[obs];
    for (int i=0;i<obs;i++) this -> G[i] = new double[obs];
    
    this -> mu = new double[obs];
    this -> sigma = new double[obs];
}

SimpleChainGMHMM::SimpleChainGMHMM(int obs, double **G, double *mu, double *sigma) : obs(obs)
{
    this -> G = new double*[obs];
    for (int i=0;i<obs;i++)
    {
        this -> G[i] = new double[obs];
        for (int j=0;j<obs;j++)
        {
            this -> G[i][j] = G[i][j];
        }
    }
    
    this -> mu = new double[obs];
    for (int i=0;i<obs;i++) this -> mu[i] = mu[i];
    
    this -> sigma = new double[obs];
    for (int i=0;i<obs;i++) this -> sigma[i] = sigma[i];
}

double gaussian_pdf(double x, double mean, double stdev)
{
    double E = x - mean;
    E *= -E;
    E /= 2 * stdev * stdev;
    double ret = exp(E);
    return ret / (stdev * sqrt(2 * M_PI));
}

double SimpleChainGMHMM::get_G(int x, int y)
{
    return G[x][y];
}

double SimpleChainGMHMM::get_probability(int obs_id, double x)
{
    return gaussian_pdf(x, mu[obs_id], sigma[obs_id]);
}

void SimpleChainGMHMM::train(vector<vector<double> > &train_set)
{
    // get means and std. deviations
    for (int i=0;i<obs;i++)
    {
        mu[i] = 0.0;
        sigma[i] = 0.0;
    }
    
    for (uint i=0;i<train_set.size();i++)
    {
        for (int j=0;j<obs;j++)
        {
            mu[j] += train_set[i][j];
        }
    }
    
    for (int i=0;i<obs;i++) mu[i] /= train_set.size();
    
    for (uint i=0;i<train_set.size();i++)
    {
        for (int j=0;j<obs;j++)
        {
            sigma[j] += (train_set[i][j] - mu[j]) * (train_set[i][j] - mu[j]);
        }
    }

    for (int i=0;i<obs;i++) sigma[i] = sqrt(sigma[i] / (train_set.size() - 1));
    
    printf("Mus/Sigmas calculated: \n");
    for (int i=0;i<obs;i++)
    {
        printf("(%lf, %lf)\n", mu[i], sigma[i]);
    }
    
    // determine observation probabilities for each position
    for (int i=0;i<obs;i++)
    {
        for (int j=0;j<obs;j++)
        {
            G[i][j] = 0.0;
        }
    }
    
    vector<vector<pair<double, int> > > sorted_vectors;
    sorted_vectors.resize(train_set.size());
    
    // initially, sort all the vectors and determine the scaling parameter
    double epsilon = 1e-6; // the smallest score possible, can be tweaked
    double max_diff = 0.0; // the maximal difference observed
    
    for (uint i=0;i<train_set.size();i++)
    {
        sorted_vectors[i].resize(obs);
        for (int j=0;j<obs;j++) sorted_vectors[i][j] = make_pair(train_set[i][j], j);
        sort(sorted_vectors[i].begin(), sorted_vectors[i].end());
        double curr_diff = sorted_vectors[i][obs-1].first - sorted_vectors[i][0].first;
        if (curr_diff > max_diff) max_diff = curr_diff;
    }
    
    double lambda = -log(epsilon) / max_diff;
    DPRINTLF(lambda);
    
    // now calculate full scores
    for (uint i=0;i<train_set.size();i++)
    {
        for (int j=0;j<obs;j++)
        {
            double pos_val = sorted_vectors[i][j].first;
            for (int k=0;k<obs;k++)
            {
                double curr_val = sorted_vectors[i][k].first;
                int curr_gene = sorted_vectors[i][k].second;
                G[j][curr_gene] += exp(-lambda * fabs(curr_val - pos_val)); // e^(-lambda*diff)
            }
        }
    }
    
    // finally, normalise the scores
    for (int i=0;i<obs;i++)
    {
        double sum = 0.0;
        for (int j=0;j<obs;j++)
        {
            sum += G[i][j];
        }
        for (int j=0;j<obs;j++)
        {
            G[i][j] /= sum;
        }
    }
    
    printf("OBSERVATION MATRIX\n");
    for (int i=0;i<obs;i++)
    {
        for (int j=0;j<obs;j++)
        {
            printf("%lf ", G[i][j]);
        }
        printf("\n");
    }
}

double SimpleChainGMHMM::log_likelihood(vector<double> &test_data)
{
    vector<pair<double, int> > sorted_data;
    sorted_data.resize(obs);
    for (int i=0;i<obs;i++) sorted_data[i] = make_pair(test_data[i], i);
    sort(sorted_data.begin(), sorted_data.end());
    
    double ret = 0.0;
    for (int i=0;i<obs;i++)
    {
        double curr_val = sorted_data[i].first;
        int curr_gene = sorted_data[i].second;
        ret += log(G[i][curr_gene]) + log(get_probability(curr_gene, curr_val));
    }
    
    return ret;
}
