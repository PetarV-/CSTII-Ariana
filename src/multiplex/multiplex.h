/*
 Petar 'PetarV' Velickovic
 Data Structure: Multiplex Network
*/

#ifndef MULTIPLEX
#define MULTIPLEX

#include <vector>

class AbstractGraphLayer
{
public:
    virtual int get_n();
    virtual double get_adj(int x, int y);
};

class SimpleGraphLayer : public AbstractGraphLayer
{
private:
    double **G;
public:
    
};

class Multiplex
{
private:
    int n;
    int L;
    std::vector<AbstractGraphLayer*> layers;
    double ****M;
public:
    Multiplex(int n, std::vector<AbstractGraphLayer*> lyrs);
    double get_edge(int lyr1, int lyr2, int n1, int n2);
    void set_edge(int lyr1, int lyr2, int n1, int n2, double val);
    double** get_matrix_form();
    double**** get_communicability_matrix();
};

#endif