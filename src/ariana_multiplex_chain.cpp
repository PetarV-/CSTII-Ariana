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

#include <hmm.h>
#include <multiplex.h>

#define KNRM  "\x1B[0m"
#define KRED  "\x1B[31m"
#define KGRN  "\x1B[32m"
#define KYEL  "\x1B[33m"
#define KBLU  "\x1B[34m"
#define KMAG  "\x1B[35m"
#define KCYN  "\x1B[36m"
#define KWHT  "\x1B[37m"

#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;
typedef unsigned int uint;
typedef long long lld;
typedef unsigned long long llu;

int main(int argc, char **argv)
{
    if (argc != 5)
    {
        printf("Usage: ./ariana_multiplex_chain <gene_count> <type_count> <training_file> <test_file>\n");
        return -1;
    }

    int gene_count = atoi(argv[1]);
    int type_count = atoi(argv[2]);
    FILE *f_train = fopen(argv[3], "r");
    FILE *f_test = fopen(argv[4], "r");
    
    printf("Importing training data...\n");
    
    int num_train;
    fscanf(f_train, "%d", &num_train);
    
    vector<vector<vector<double> > > train_patient, train_normal;
    int patient_cnt = 0, normal_cnt = 0;
    
    for (int i=0;i<num_train;i++)
    {
        char status[10];
        fscanf(f_train, "%s", status);
        if (strcmp(status, "patient") == 0)
        {
            train_patient.push_back(vector<vector<double> >(gene_count));
            train_patient[patient_cnt].resize(gene_count);
            for (int j=0;j<gene_count;j++) train_patient[patient_cnt][j].resize(type_count);
            for (int k=0;k<type_count;k++)
            {
                for (int j=0;j<gene_count;j++)
                {
                    fscanf(f_train, "%lf", &train_patient[patient_cnt][j][k]);
                }
            }
            patient_cnt++;
        }
        else if (strcmp(status, "normal") == 0)
        {
            train_normal.push_back(vector<vector<double> >(gene_count));
            train_normal[normal_cnt].resize(gene_count);
            for (int j=0;j<gene_count;j++) train_normal[normal_cnt][j].resize(type_count);
            for (int k=0;k<type_count;k++)
            {
                for (int j=0;j<gene_count;j++)
                {
                    fscanf(f_train, "%lf", &train_normal[normal_cnt][j][k]);
                }
            }
            normal_cnt++;
        }
    }
    
    printf("Training data imported!\n");
    
    printf("Training the models...\n");
    
    HMMChainMultiplex *patient_model = new HMMChainMultiplex(gene_count, type_count);
    HMMChainMultiplex *normal_model = new HMMChainMultiplex(gene_count, type_count);
    
    patient_model -> train(train_patient);
    normal_model -> train(train_normal);
    
    printf("Models trained!\n");
    
    printf("Classifying test data...\n");
    
    int num_test;
    fscanf(f_test, "%d", &num_test);
    
    int true_patients = 0, false_patients = 0;
    int true_normal = 0, false_normal = 0;
    
    for (int i=0;i<num_test;i++)
    {
        printf("--------------------------------------------\n");
        printf("Classifying test vector #%d... ", i);
        
        char status[10];
        bool expected_inference;
        vector<vector<double> > test_vector;
        test_vector.resize(gene_count);
        for (int j=0;j<gene_count;j++) test_vector[j].resize(type_count);
        
        fscanf(f_test, "%s", status);
        if (strcmp(status, "patient") == 0) expected_inference = true;
        else if (strcmp(status, "normal") == 0) expected_inference = false;
        
        for (int k=0;k<type_count;k++)
        {
            for (int j=0;j<gene_count;j++)
            {
                fscanf(f_test, "%lf", &test_vector[j][k]);
            }
        }
        
        double lhood1 = patient_model -> log_likelihood(test_vector);
        double lhood0 = normal_model -> log_likelihood(test_vector);
        
        printf("OK!\n");
        
        printf("likelihood(normal) = %lf, likelihood(patient) = %lf\n", lhood0, lhood1);
        
        bool inference = (lhood1 > lhood0);
        
        if (inference && expected_inference) true_patients++;
        else if (!inference && expected_inference) false_normal++;
        else if (inference && !expected_inference) false_patients++;
        else if (!inference && !expected_inference) true_normal++;
        
        printf("Test vector classified as %s... %s\n", inference ? "patient" : "normal", (inference == expected_inference) ? (KGRN "CORRECT!" KNRM) : (KRED "INCORRECT!" KNRM));
    }
    
    printf("--------------------------------------------\n");
    
    printf("============================================\n");
    printf("Classifier evaluation results:\n");
    printf("Precision: %lf%%\n", 100.0 * true_patients / (1.0 * (true_patients + false_patients)));
    printf("Sensitivity: %lf%%\n", 100.0 * true_patients / (1.0 * (true_patients + false_normal)));
    printf("Specificity: %lf%%\n", 100.0 * true_normal / (1.0 * (true_normal + false_patients)));
    printf("Accuracy: %lf%%\n", 100.0 * (true_patients + true_normal) / (1.0 * num_test));
    printf("============================================\n");
    
    fclose(f_train);
    fclose(f_test);
    
    printf("Done.\n");
}