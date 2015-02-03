#ifndef CLASSIFIER
#define CLASSIFIER

#include <vector>

#include <hmm.h>
#include <multiplex.h>

template<typename Data, typename Label>
class Classifier
{
public:
    virtual ~Classifier() { }
    virtual void train(std::vector<std::pair<Data, Label> > &training_set) = 0;
    virtual Label classify(Data &test_data) = 0;
    virtual std::vector<std::pair<double, Label> > get_thresholds() = 0;
};

class SingleChainClassifier : public Classifier<std::vector<std::vector<double> >, bool>
{
private:
    int gene_count;
    int param_id;
    SimpleChainGMHMM* patient_model;
    SimpleChainGMHMM* normal_model;
    
    std::vector<std::pair<double, bool> > thresholds;
    
public:
    SingleChainClassifier(int gene_count, int param_id = 0);
    ~SingleChainClassifier();
    
    void train(std::vector<std::pair<std::vector<std::vector<double> >, bool> > &training_set);
    bool classify(std::vector<std::vector<double> > &test_data);
    
    std::vector<std::pair<double, bool> > get_thresholds();
};

class MultiplexChainClassifier : public Classifier<std::vector<std::vector<double> >, bool>
{
private:
    int gene_count;
    int type_count;
    HMMChainMultiplex* patient_model;
    HMMChainMultiplex* normal_model;
    
    std::vector<std::pair<double, bool> > thresholds;
    
public:
    MultiplexChainClassifier(int gene_count, int type_count);
    ~MultiplexChainClassifier();
    
    void train(std::vector<std::pair<std::vector<std::vector<double> >, bool> > &training_set);
    bool classify(std::vector<std::vector<double> > &test_data);
    
    std::vector<std::pair<double, bool> > get_thresholds();
};

#endif
