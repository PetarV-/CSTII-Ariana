/*
 Petar 'PetarV' Velickovic
 Data Structure: Hidden Markov Model
*/

#ifndef HIDDEN_MARKOV_MODEL
#define HIDDEN_MARKOV_MODEL

#include <vector>

class HMM
{
private:
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

    void forward_backward(std::vector<int> Y);
    std::vector<int> viterbi(std::vector<int> Y);
    void baumwelch(std::vector<int> Y);
};

#endif
