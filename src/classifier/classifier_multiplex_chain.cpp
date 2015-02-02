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

#include <classifier.h>
#include <multiplex.h>

#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;
typedef unsigned int uint;
typedef long long lld;
typedef unsigned long long llu;

MultiplexChainClassifier::MultiplexChainClassifier(int gene_count, int type_count) : gene_count(gene_count), type_count(type_count)
{
    patient_model = new HMMChainMultiplex(gene_count, type_count);
    normal_model = new HMMChainMultiplex(gene_count, type_count);
}

MultiplexChainClassifier::~MultiplexChainClassifier()
{
    delete patient_model;
    delete normal_model;
}

void MultiplexChainClassifier::train(vector<pair<vector<vector<double> >, bool> > &training_set)
{
    vector<vector<vector<double> > > train_patient, train_normal;
    int patient_cnt = 0, normal_cnt = 0;
    
    for (uint i=0;i<training_set.size();i++)
    {
        if (training_set[i].second)
        {
            train_patient.push_back(vector<vector<double> >(gene_count));
            train_patient[patient_cnt].resize(gene_count);
            for (int j=0;j<gene_count;j++)
            {
                train_patient[patient_cnt][j].resize(type_count);
                for (int k=0;k<type_count;k++)
                {
                    train_patient[patient_cnt][j][k] = training_set[i].first[j][k];
                }
            }
            patient_cnt++;
        }
        else
        {
            train_normal.push_back(vector<vector<double> >(gene_count));
            train_normal[normal_cnt].resize(gene_count);
            for (int j=0;j<gene_count;j++)
            {
                train_normal[normal_cnt][j].resize(type_count);
                for (int k=0;k<type_count;k++)
                {
                    train_patient[normal_cnt][j][k] = training_set[i].first[j][k];
                }
            }
            normal_cnt++;
        }
    }
    
    patient_model -> train(train_patient);
    normal_model -> train(train_normal);
}

bool MultiplexChainClassifier::classify(vector<vector<double> > &test_data)
{
    double lhood1 = patient_model -> log_likelihood(test_data);
    double lhood0 = normal_model -> log_likelihood(test_data);
    
    return (lhood1 > lhood0);
}
