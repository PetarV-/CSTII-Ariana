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
    ~HMM();
    
    void forward_backward(std::vector<int> &Y);
    std::vector<int> viterbi(std::vector<int> &Y);
    void baumwelch(std::vector<int> &Y);
};

class SimpleChainGMHMM
{
private:
    int obs; // number of states and "sub-observations"
    double **G; // sub-observation emission probabilites
    double *mu, *sigma; // means and variances for each sub-observation
    
public:
    SimpleChainGMHMM(int obs);
    SimpleChainGMHMM(int obs, double **G, double *mu, double *sigma);
    ~SimpleChainGMHMM();
    
    double get_G(int x, int y);
    double get_probability(int obs_id, double x);
    void train(std::vector<std::vector<double> > &train_set);
    double log_likelihood(std::vector<double> &test_data);
};

#endif
