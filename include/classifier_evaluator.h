#ifndef CLASSIFIER_EVAL
#define CLASSIFIER_EVAL

#include <vector>

#include <classifier.h>

struct run_result
{
    int true_positives, false_positives;
    int false_negatives, true_negatives;
    
    double accuracy;
    double precision;
    double sensitivity;
    double specificity;
    double false_positive_rate;
    double negative_predictive_value;
    double false_discovery_rate;
    
    double mcc;
    double f1_score;
    
    std::vector<std::pair<double, bool> > likelihood_diffs;
};

template<typename Data, typename Label>
class Classifier
{
public:
    virtual ~Classifier() { }
    virtual void train(std::vector<std::pair<Data, Label> > &training_set) = 0;
    virtual Label classify(Data &test_data) = 0;
};

class SingleChainClassifier : public Classifier<std::vector<std::vector<double> >, bool>
{
private:
    int gene_count;
    int param_id;
    SimpleChainGMHMM* patient_model;
    SimpleChainGMHMM* normal_model;
    
public:
    SingleChainClassifier(int gene_count, int param_id = 0);
    ~SingleChainClassifier();
    
    void train(std::vector<std::pair<std::vector<std::vector<double> >, bool> > &training_set);
    bool classify(std::vector<std::vector<double> > &test_data);
};

class MultiplexChainClassifier : public Classifier<std::vector<std::vector<double> >, bool>
{
private:
    int gene_count;
    int type_count;
    HMMChainMultiplex* patient_model;
    HMMChainMultiplex* normal_model;
    
public:
    MultiplexChainClassifier(int gene_count, int type_count);
    ~MultiplexChainClassifier();
    
    void train(std::vector<std::pair<std::vector<std::vector<double> >, bool> > &training_set);
    bool classify(std::vector<std::vector<double> > &test_data);
};

#endif
