/*
 Petar 'PetarV' Velickovic
 Data Structure: Abstract Graph Layer (+ Implementations)
*/

#ifndef LAYERS
#define LAYERS

class AbstractGraphLayer
{
public:
    virtual int get_n() = 0;
    virtual double get_adj(int x, int y) = 0;
    virtual void set_adj(int x, int y, double val) = 0;
};

class SimpleGraphLayer : public AbstractGraphLayer
{
private:
    int n;
    double **G;
public:
    SimpleGraphLayer(int n);
    SimpleGraphLayer(int n, double **A);
    int get_n();
    double get_adj(int x, int y);
    void set_adj(int x, int y, double val);
};

#endif