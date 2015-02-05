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
    if (argc != 3)
    {
        printf("Usage: ./simple_tester <classifier_type> <data_set_file>\n");
        return -1;
    }
    
    string classifier_type = argv[1];
    vector<pair<vector<vector<double> >, bool> > data = extract_data(argv[2]);
    
    int gene_count = data[0].first.size();
    int type_count = data[0].first[0].size();
    
    if (classifier_type == "gmhmm-single-chain") C = new SingleChainClassifier(gene_count);
    else if (classifier_type == "gmhmm-multiplex-chain") C = new MultiplexChainClassifier(gene_count, type_count);
    else if (classifier_type == "gmhmm-single");
    else if (classifier_type == "gmhmm-multiplex");
    else if (classifier_type == "general-multiplex");
    else
    {
        printf("classifier_type must be of the following kind:\n");
        printf("gmhmm-single-chain\n");
        printf("gmhmm-multiplex-chain\n");
        printf("gmhmm-single\n");
        printf("gmhmm-multiplex\n");
        printf("general-multiplex\n");
        return -2;
    }
    
    run_result results = crossvalidate(C, data);
    
    return 0;
}