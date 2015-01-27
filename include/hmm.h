/*
 Petar 'PetarV' Velickovic
 Data Structure: Hidden Markov Model
*/

#ifndef HIDDEN_MARKOV_MODEL
#define HIDDEN_MARKOV_MODEL

#include <vector>

class HMM
{
protected:
    int n; // number of nodes
    int obs; // number of observations
    double **T; // transition probability matrix
    double **P; // output probability matrix
    
    double **A, **B, **Pst; // alpha, beta, state prob. matrices
    int *Arn, *Brn; // how many times have we renormalised A, B
    int renorm;
    double likelihood;

public:
    HMM(int n, int obs);
    HMM(int n, int obs, double **T, double **P);

    void forward_backward(std::vector<int> &Y);
    std::vector<int> viterbi(std::vector<int> &Y);
    void baumwelch(std::vector<int> &Y);
};

class SimpleChainGMHMM
{
private:
    int obs; // number of observations (also nodes)
    double **G; // observation probabilites
    double *mu, *sigma; // means and variances of each observation
    
public:
    SimpleChainGMHMM(int obs);
    SimpleChainGMHMM(int obs, double **G, double *mu, double *sigma);
    
    double get_G(int x, int y);
    double get_probability(int obs_id, double x);
    void train(std::vector<std::vector<double> > &train_set);
    double log_likelihood(std::vector<double> &test_data);
};

class tdHMM : public HMM
{
private:
    int k; // the dimension of the HMM
    double *P_init; // initial probability of each dimension (k)
    double **dim_trans; // transition probabilities between dimensions (k x k)
    double **dim_obs; // observation probabilities of each dimension (k x obs)
public:
    tdHMM(int n, int k, int obs, double *P, double **A, double **B);
};

#endif
