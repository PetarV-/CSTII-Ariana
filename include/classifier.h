#ifndef CLASSIFIER
#define CLASSIFIER

#include <vector>

#include <hmm.h>
#include <layers.h>
#include <multiplex.h>

template<typename Data, typename Label>
class Classifier
{
public:
    virtual ~Classifier() { }
    virtual void train(std::vector<std::pair<Data, Label> > &training_set) = 0;
    virtual Label classify(Data &test_data) = 0;
    virtual std::vector<double> get_thresholds() = 0;
};

class SingleChainClassifier : public Classifier<std::vector<std::vector<double> >, bool>
{
private:
    int gene_count;
    int param_id;
    SimpleChainGMHMM* patient_model;
    SimpleChainGMHMM* normal_model;
    
    std::vector<double> thresholds;
    
public:
    SingleChainClassifier(int gene_count, int param_id = 0);
    ~SingleChainClassifier();
    
    void train(std::vector<std::pair<std::vector<std::vector<double> >, bool> > &training_set);
    bool classify(std::vector<std::vector<double> > &test_data);
    
    std::vector<double> get_thresholds();
};

class MultiplexChainClassifier : public Classifier<std::vector<std::vector<double> >, bool>
{
private:
    int gene_count;
    int type_count;
    HMMChainMultiplex* patient_model;
    HMMChainMultiplex* normal_model;
    
    std::vector<double> thresholds;
    
public:
    MultiplexChainClassifier(int gene_count, int type_count);
    ~MultiplexChainClassifier();
    
    void train(std::vector<std::pair<std::vector<std::vector<double> >, bool> > &training_set);
    bool classify(std::vector<std::vector<double> > &test_data);
    
    std::vector<double> get_thresholds();
};

class GMHMMClassifier : public Classifier<std::vector<std::vector<double> >, bool>
{
private:
    int gene_count;
    int param_id;
    GMHMM* patient_model;
    GMHMM* normal_model;
    
    std::vector<double> thresholds;
    
public:
    GMHMMClassifier(int gene_count, int param_id = 0);
    ~GMHMMClassifier();
    
    void train(std::vector<std::pair<std::vector<std::vector<double> >, bool> > &training_set);
    bool classify(std::vector<std::vector<double> > &test_data);
    
    std::vector<double> get_thresholds();
};

class MultiplexGMHMMClassifier : public Classifier<std::vector<std::vector<double> >, bool>
{
private:
    int gene_count;
    int type_count;
    MultiplexGMHMM* patient_model;
    MultiplexGMHMM* normal_model;
    
    std::vector<double> thresholds;
    
public:
    MultiplexGMHMMClassifier(int gene_count, int type_count);
    ~MultiplexGMHMMClassifier();
    
    void train(std::vector<std::pair<std::vector<std::vector<double> >, bool> > &training_set);
    bool classify(std::vector<std::vector<double> > &test_data);
    
    std::vector<double> get_thresholds();
};

class GenericSingleLayerClassifier : public Classifier<std::vector<std::vector<double> >, bool>
{
private:
    int gene_count;
    int param_id;
    SimpleGraphLayer* patient_model;
    SimpleGraphLayer* normal_model;
    
    std::vector<double> thresholds;
    
public:
    GenericSingleLayerClassifier(int gene_count, int param_id = 0);
    ~GenericSingleLayerClassifier();
    
    void train(std::vector<std::pair<std::vector<std::vector<double> >, bool> > &training_set);
    bool classify(std::vector<std::vector<double> > &test_data);
    
    std::vector<double> get_thresholds();
};

class GenericMultiplexClassifier : public Classifier<std::vector<std::vector<double> >, bool>
{
private:
    int gene_count;
    int type_count;
    Multiplex* patient_model;
    Multiplex* normal_model;
    
    std::vector<double> thresholds;
    
public:
    GenericMultiplexClassifier(int gene_count, int type_count);
    ~GenericMultiplexClassifier();
    
    void train(std::vector<std::pair<std::vector<std::vector<double> >, bool> > &training_set);
    bool classify(std::vector<std::vector<double> > &test_data);
    
    std::vector<double> get_thresholds();
};

#endif
