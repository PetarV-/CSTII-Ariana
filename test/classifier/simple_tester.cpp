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

Classifier<vector<vector<double> >, bool> *C;

int main(int argc, char **argv)
{
    if (argc != 5)
    {
        printf("Usage: ./simple_tester <classifier_type> <data_set_file> <noise_mean_lo>:<noise_mean_step>:<noise_mean_hi> <noise_stddev_lo>:<noise_stddev_step>:<noise_stddev_hi>\n");
        return -1;
    }
    
    string classifier_type = argv[1];
    vector<pair<vector<vector<double> >, bool> > data = extract_data(argv[2]);
    double mu_lo, mu_step, mu_hi;
    double sigma_lo, sigma_step, sigma_hi;
    sscanf(argv[3], "%lf:%lf:%lf", &mu_lo, &mu_step, &mu_hi);
    sscanf(argv[4], "%lf:%lf:%lf", &sigma_lo, &sigma_step, &sigma_hi);
    
    int gene_count = data[0].first.size();
    int type_count = data[0].first[0].size();
    
    if (classifier_type == "gmhmm-single-chain") C = new SingleChainClassifier(gene_count);
    else if (classifier_type == "gmhmm-multiplex-chain") C = new MultiplexChainClassifier(gene_count, type_count);
    else if (classifier_type == "gmhmm-single") C = new GMHMMClassifier(gene_count);
    else if (classifier_type == "gmhmm-multiplex") C = new MultiplexGMHMMClassifier(gene_count, type_count);
    else if (classifier_type == "general-single") C = new GenericSingleLayerClassifier(gene_count);
    else if (classifier_type == "general-multiplex") C = new GenericMultiplexClassifier(gene_count, type_count);
    else
    {
        printf("classifier_type must be of the following kind:\n");
        printf("gmhmm-single-chain\n");
        printf("gmhmm-multiplex-chain\n");
        printf("gmhmm-single\n");
        printf("gmhmm-multiplex\n");
        printf("general-single\n");
        printf("general-multiplex\n");
        return -2;
    }
    
    noise_test(C, data, mu_lo, mu_step, mu_hi, sigma_lo, sigma_step, sigma_hi);
    
    return 0;
}
