/*
 Petar 'PetarV' Velickovic
 Data Structure: Abstract Graph Layer (+ Sample Implementation)
*/

#ifndef LAYERS
#define LAYERS

class AbstractGraphLayer
{
public:
    virtual ~AbstractGraphLayer() { }
    virtual int get_n() = 0;
    virtual double get_adj(int x, int y) = 0;
    virtual void set_adj(int x, int y, double val) = 0;
    
    virtual void train(std::vector<std::vector<double> > &train_set) = 0;
    virtual double log_likelihood(std::vector<double> &test_data) = 0;
};

class SimpleGraphLayer : public AbstractGraphLayer
{
private:
    int n;
    double **G;
public:
    SimpleGraphLayer(int n);
    SimpleGraphLayer(int n, double **A);
    ~SimpleGraphLayer();
    
    int get_n();
    double get_adj(int x, int y);
    void set_adj(int x, int y, double val);
    
    void train(std::vector<std::vector<double> > &train_set);
    double log_likelihood(std::vector<double> &test_data);
};

#endif
