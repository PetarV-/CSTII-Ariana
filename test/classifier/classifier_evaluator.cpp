#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <deque>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <complex>
#include <functional>

#include <classifier.h>
#include <classifier_evaluator.h>

#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;
typedef unsigned int uint;
typedef long long lld;
typedef unsigned long long llu;

run_result single_run(Classifier<vector<vector<double> >, bool> *C, vector<pair<vector<vector<double> >, bool> > &training_set, vector<pair<vector<vector<double> >, bool> > &test_set)
{
    run_result ret;
    ret.true_positives = ret.false_positives = 0;
    ret.false_negatives = ret.true_negatives = 0;
    
    C -> train(training_set);
    
    int total = test_set.size();
    int total_positives = 0, total_negatives = 0;
    
    for (uint i=0;i<test_set.size();i++)
    {
        bool expected_inference = test_set[i].second;
        bool inference = C -> classify(test_set[i].first);
        
        if (inference && expected_inference) ret.true_positives++;
        else if (inference && !expected_inference) ret.false_positives++;
        else if (!inference && expected_inference) ret.false_negatives++;
        else if (!inference && !expected_inference) ret.true_negatives++;
        
        if (expected_inference) total_positives++;
        else total_negatives++;
    }
    
    ret.accuracy = (ret.true_positives + ret.true_negatives) * 1.0 / total * 1.0;
    ret.precision = ret.true_positives * 1.0 / (ret.true_positives + ret.false_positives) * 1.0;
    ret.sensitivity = ret.true_positives * 1.0 / (ret.true_positives + ret.false_negatives) * 1.0;
    ret.specificity = ret.true_negatives * 1.0 / (ret.true_negatives + ret.false_positives) * 1.0;
    ret.false_positive_rate = ret.false_positives * 1.0 / (ret.false_positives + ret.true_negatives) * 1.0;
    ret.negative_predictive_value = ret.true_negatives * 1.0 / (ret.true_negatives + ret.false_negatives) * 1.0;
    ret.false_discovery_rate = ret.false_positives * 1.0 / (ret.false_positives + ret.true_positives) * 1.0;
    
    ret.mcc = (ret.true_positives * ret.true_negatives - ret.false_positives * ret.false_negatives) * 1.0 / sqrt((ret.true_positives + ret.false_positives) * (ret.true_negatives + ret.false_negatives) * (ret.true_positives + ret.false_negatives) * (ret.true_negatives + ret.false_positives));
    ret.f1_score = ret.precision * ret.sensitivity * 2.0 / (ret.precision + ret.sensitivity) * 1.0;
    
    vector<double> thresh = C -> get_thresholds();
    vector<pair<double, bool> > roc_meta;
    roc_meta.resize(test_set.size());
    for (uint i=0;i<test_set.size();i++)
    {
        roc_meta[i] = make_pair(thresh[i], test_set[i].second);
    }
    sort(roc_meta.begin(), roc_meta.end(), greater<pair<double, bool> >());
    
    ret.roc_points.resize(total + 1);
    ret.roc_points[0] = make_pair(roc_meta[0].first + 1.0, make_pair(0.0, 0.0));
    ret.roc_auc = 0.0;
    
    int curr_true_positives = 0, curr_false_positives = 0;
    int curr_false_negatives = total_positives, curr_true_negatives = total_negatives;
    
    for (uint i=0;i<roc_meta.size();i++)
    {
        if (roc_meta[i].second) curr_true_positives++, curr_false_negatives--;
        else curr_false_positives++, curr_true_negatives--;
        
        double old_sensitivity = ret.roc_points[i].second.first;
        double old_fpr = ret.roc_points[i].second.second;
        
        double new_sensitivity = curr_true_positives * 1.0 / (curr_true_positives + curr_false_negatives) * 1.0;
        double new_fpr = curr_false_positives * 1.0 / (curr_false_positives + curr_true_negatives) * 1.0;
        
        if (!roc_meta[i].second) ret.roc_auc += old_sensitivity * (new_fpr - old_fpr);
        ret.roc_points[i+1] = make_pair(roc_meta[i].first, make_pair(new_sensitivity, new_fpr));
    }
    
    return ret;
}

run_result crossvalidate(Classifier<vector<vector<double> >, bool> *C, vector<pair<vector<vector<double> >, bool> > &training_set, int fold_cnt)
{
    int total = training_set.size();
    int total_patient = 0, total_normal = 0;
    for (uint i=0;i<training_set.size();i++)
    {
        if (training_set[i].second) total_patient++;
        else total_normal++;
    }
    
    int fold_size_patient = total_patient / fold_cnt;
    int fold_size_normal = total_normal / fold_cnt;
    int rem_patient = total_patient % fold_cnt;
    int rem_normal = total_normal % fold_cnt;
    
    vector<vector<pair<vector<vector<double> >, bool> > > folds;
    folds.resize(fold_cnt);
    
    int *fold_size = new int[fold_cnt];
    for (int i=0;i<fold_cnt;i++)
    {
        fold_size[i] = fold_size_patient + (i < rem_patient) + fold_size_normal + (i < rem_normal);
        folds[i].resize(fold_size[i]);
    }
    
    int curr_patient_fold = 0, curr_normal_fold = 0;
    int curr_patient_fold_size = fold_size_patient + (rem_patient > 0);
    int curr_normal_fold_size = fold_size_normal + (rem_normal > 0);
    
    for (uint i=0;i<training_set.size();i++)
    {
        if (training_set[i].second)
        {
            folds[curr_patient_fold][--fold_size[curr_patient_fold]].second = training_set[i].second;
            folds[curr_patient_fold][fold_size[curr_patient_fold]].first.resize(training_set[i].first.size());
            copy(training_set[i].first.begin(), training_set[i].first.end(), folds[curr_patient_fold][fold_size[curr_patient_fold]].first.begin());
            if (--curr_patient_fold_size == 0) curr_patient_fold_size = fold_size_patient + (++curr_patient_fold < rem_patient);
        }
        else
        {
            folds[curr_normal_fold][--fold_size[curr_normal_fold]].second = training_set[i].second;
            folds[curr_normal_fold][fold_size[curr_normal_fold]].first.resize(training_set[i].first.size());
            copy(training_set[i].first.begin(), training_set[i].first.end(), folds[curr_normal_fold][fold_size[curr_normal_fold]].first.begin());
            if (--curr_normal_fold_size == 0) curr_normal_fold_size = fold_size_normal + (++curr_normal_fold < rem_normal);
        }
    }
    
    delete[] fold_size;
    
    vector<run_result> individual;
    individual.resize(fold_cnt);
    for (int i=0;i<fold_cnt;i++)
    {
        vector<pair<vector<vector<double> >, bool> > curr_train, curr_test;
        curr_test = folds[i];
        curr_train.resize(total - folds[i].size());
        for (int j=0;j<fold_cnt;j++)
        {
            if (i != j) curr_train.insert(curr_train.end(), folds[j].begin(), folds[j].end());
        }
        
        individual[i] = single_run(C, curr_train, curr_test);
    }

    // TODO
    return individual[0];
}

vector<pair<vector<vector<double> >, bool> > extract_data(char* filename)
{
    int total;
    int gene_count, type_count;
    char expected_outcome[101];
    
    vector<pair<vector<vector<double> >, bool> > ret;
    
    FILE *f = fopen(filename, "r");
    
    fscanf(f, "%d", &total);
    fscanf(f, "%d%d", &gene_count, &type_count);
    
    ret.resize(total);
    
    for (int i=0;i<total;i++)
    {
        ret[i].first.resize(gene_count);
        for (int j=0;j<gene_count;j++) ret[i].first[j].resize(type_count);
        
        fscanf(f, "%s", expected_outcome);
        for (int k=0;k<type_count;k++)
        {
            for (int j=0;j<gene_count;j++)
            {
                fscanf(f, "%lf", &ret[i].first[j][k]);
            }
        }
        ret[i].second = (strcmp(expected_outcome, "patient") == 0);
    }
    
    fclose(f);
    
    return ret;
}
