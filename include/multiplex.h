/*
 Petar 'PetarV' Velickovic
 Data Structure: Multiplex Network
*/

#ifndef MULTIPLEX
#define MULTIPLEX

#include <functional>
#include <vector>

#include <layers.h>
#include <hmm.h>

class Multiplex
{
private:
    int n;
    int L;
    std::vector<AbstractGraphLayer*> layers;
    double ****M;
    
    std::vector<std::function<double(std::vector<double>)> > objectives;
public:
    Multiplex(int n, int L);
    Multiplex(int n, int L, double ****A);
    Multiplex(int n, std::vector<AbstractGraphLayer*> &lyrs);
    Multiplex(int n, std::vector<AbstractGraphLayer*> &lyrs, double **omega);
    Multiplex(int n, std::vector<AbstractGraphLayer*> &lyrs, double ****inter_lyr);
    ~Multiplex();
    
    double get_edge(int lyr1, int lyr2, int n1, int n2);
    void set_edge(int lyr1, int lyr2, int n1, int n2, double val);
    double** get_matrix_form();
    double**** get_communicability_matrix();
    double** get_aggregate_matrix(bool normalise);
    
    void sync();
    void set_omega(double **omega, bool normalise, bool upsync);
    void train(std::vector<std::vector<std::vector<double> > > &train_set);
    double log_likelihood(std::vector<std::vector<double> > &test_data);
    
    void dump_muxviz_data(char *nodes_filename, char *base_layers_filename);
};

class HMMChainMultiplex
{
private:
    int obs;
    int L;
    std::vector<SimpleChainGMHMM*> layers;
    double **omega;
    
    std::vector<std::function<double(std::vector<double>)> > objectives;
public:
    HMMChainMultiplex(int obs, int L);
    ~HMMChainMultiplex();
    
    void set_omega(double **omega);
    void train(std::vector<std::vector<std::vector<double> > > &train_set);
    double log_likelihood(std::vector<std::vector<double> > &test_data);
    
    void dump_muxviz_data(char *nodes_filename, char *base_layers_filename);
};

class MultiplexGMHMM
{
private:
    int n, obs;
    int L;
    std::vector<GMHMM*> layers;
    double **omega;
    
    std::vector<std::function<double(std::vector<double>)> > objectives;
public:
    MultiplexGMHMM(int n, int obs, int L);
    ~MultiplexGMHMM();
    
    void set_omega(double **omega);
    void train(std::vector<std::vector<std::vector<double> > > &train_set);
    double log_likelihood(std::vector<std::vector<double> > &test_data);
    
    void dump_muxviz_data(char *nodes_filename, char *base_layers_filename);
};

#endif
