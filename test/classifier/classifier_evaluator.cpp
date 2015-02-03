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

#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;
typedef unsigned int uint;
typedef long long lld;
typedef unsigned long long llu;

struct single_run_results

int main(int argc, char **argv)
{
    if (argc < 4 || argc > 5)
    {
        printf("Usage: ./classifier_evaluator <classifier_type> <gene_count> <type_count> <data_set_file>\n");
        return -1;
    }
    
    int gene_count = atoi(argv[1]);
    int type_count = atoi(argv[2]);
    FILE *f_train = fopen(argv[3], "r");
    FILE *f_test = fopen(argv[4], "r");
    
    printf("Importing training data...\n");
    
    int num_train;
    fscanf(f_train, "%d", &num_train);
}
