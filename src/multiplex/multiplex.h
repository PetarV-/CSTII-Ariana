/*
 Petar 'PetarV' Velickovic
 Data Structure: Multiplex Network
*/

#ifndef MULTIPLEX
#define MULTIPLEX

#include <vector>

#include "../layers/layers.h"
#include "../hmm/hmm.h"

class Multiplex
{
private:
    int n;
    int L;
    std::vector<AbstractGraphLayer*> layers;
    double ****M;
public:
    Multiplex(int n, int L);
    Multiplex(int n, int L, double ****A);
    Multiplex(int n, std::vector<AbstractGraphLayer*> lyrs);
    Multiplex(int n, std::vector<AbstractGraphLayer*> lyrs, double **omega);
    Multiplex(int n, std::vector<AbstractGraphLayer*> lyrs, double ****inter_lyr);
    double get_edge(int lyr1, int lyr2, int n1, int n2);
    void set_edge(int lyr1, int lyr2, int n1, int n2, double val);
    double** get_matrix_form();
    double**** get_communicability_matrix();
    double** get_aggregate_matrix();
};

class HMMChainMultiplex
{
private:
    int obs;
    int L;
    std::vector<SimpleChainGMHMM*> layers;
    double **omega;
public:
    HMMChainMultiplex(int obs, int L);
    
    void set_omega(double **omega);
    void train(std::vector<std::vector<std::vector<double> > > train_set);
    double log_likelihood(std::vector<std::vector<double> > test_data);
};

#endif
