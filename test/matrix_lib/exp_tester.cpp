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

#include "../../src/matrix_lib/matrix_lib.h"

#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;
typedef long long lld;
typedef unsigned long long llu;

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        printf("Usage: ./exp_tester <input_matrix_file> <output_file>\n");
        return -1;
    }
    
    int len;
    double **mat;
    
    FILE *f = fopen(argv[1], "r");
    
    fscanf(f, "%d", &len);
    
    mat = new double*[len];
    
    for (int i=0;i<len;i++)
    {
        mat[i] = new double[len];
        for (int j=0;j<len;j++)
        {
            fscanf(f, "%lf", &mat[i][j]);
        }
    }
    
    FILE *g = fopen(argv[2], "w");
    
    fprintf(g, "EXP(M) = \n");
    double **fin = mat_exp(mat, len);
    for (int i=0;i<len;i++)
    {
        for (int j=0;j<len;j++)
        {
            fprintf(g, "%lf ", fin[i][j]);
        }
        fprintf(g, "\n");
    }
    
    fclose(f);
    fclose(g);
    
    printf("Done. Matrix exponential written to %s.\n", argv[2]);
    
    return 0;
}
