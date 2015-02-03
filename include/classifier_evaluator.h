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
    
    std::vector<std::pair<double, double> > roc_points;
    double roc_auc;
};

run_result single_run(Classifier<std::vector<std::vector<double> >, bool> *C, std::vector<std::pair<std::vector<std::vector<double> >, bool> > &training_set, std::vector<std::pair<std::vector<std::vector<double> >, bool> > &test_set);

run_result crossvalidate(Classifier<std::vector<std::vector<double> >, bool> *C, std::vector<std::pair<std::vector<std::vector<double> >, bool> > &training_set, int folds = 10);

std::vector<std::pair<std::vector<std::vector<double> >, bool> > extract_data(char* filename);

#endif
