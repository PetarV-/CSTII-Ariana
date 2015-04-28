/*
 Petar 'PetarV' Velickovic
 Data Structure: Hidden Markov Model
*/

#ifndef HIDDEN_MARKOV_MODEL
#define HIDDEN_MARKOV_MODEL

#include <tuple>
#include <vector>

// helper function; Phi(x; mean, stddev)
double gaussian_pdf(double x, double mean, double stdev);

class HMM
{
protected:
    int n; // number of nodes
    int obs; // number of observations
    double *pi; // start-state probability vector
    double **T; // transition probability matrix
    double **P; // output probability matrix
    
    double **A, **B, **Pst; // alpha, beta, state prob. matrices
    int *Arn, *Brn; // how many times have we renormalised A, B
    int renorm;
    double likelihood;

public:
    HMM(int n, int obs);
    HMM(int n, int obs, double *pi, double **T, double **P);
    ~HMM();
    
    double** get_A();
    double** get_B();
    double* get_pi();
    double** get_T();
    double** get_P();
    
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

class GMHMM
{
private:
    int n, obs; // number of nodes and "sub-observations"
    double *pi; // start-state probability vector
    double **T; // transition probability matrix
    double **O; // sub-output emission matrix
    double *mu, *sigma; // means and variances for each sub-output
    
public:
    GMHMM(int n, int obs); // initialise a random GMHMM
    GMHMM(int n, int obs, double *pi, double **T, double **O, double *mu, double *sigma); // load a known GMHMM
    ~GMHMM();
    
    std::tuple<double**, double*, double> forward(std::vector<std::pair<double, int> > &Y);
    double** backward(std::vector<std::pair<double, int> > &Y, double *c);
    void baumwelch(std::vector<std::vector<double> > &Ys, int iterations, double tolerance);
    
    double get_pi(int x);
    double get_T(int i, int j);
    double get_O(int x, int y);
    double get_probability(int obs_id, double x);
    void train(std::vector<std::vector<double> > &train_set);
    double log_likelihood(std::vector<double> &test_data);
};

#endif
