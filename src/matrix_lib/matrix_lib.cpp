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

#include "matrix_lib.h"

#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;
typedef long long lld;
typedef unsigned long long llu;

void solve(double **A, double **B, int n, int nb)
{
    for (int cur_j=0;cur_j<n;cur_j++)
    {
        double piv = fabs(A[cur_j][cur_j]);
        int best = cur_j;
        for (int i=cur_j+1;i<n;i++)
        {
            if (fabs(A[i][cur_j]) > piv)
            {
                piv = fabs(A[i][cur_j]);
                best = i;
            }
        }
        
        if (best != cur_j)
        {
            for (int j=0;j<n;j++) swap(A[cur_j][j], A[best][j]);
            for (int j=0;j<nb;j++) swap(B[cur_j][j], B[best][j]);
        }
        
        double t = A[cur_j][cur_j];
        A[cur_j][cur_j] = 1.0;
        for (int j=cur_j+1;j<n;j++) A[cur_j][j] /= t;
        for (int j=0;j<nb;j++) B[cur_j][j] /= t;
        
        for (int i=cur_j+1;i<n;i++)
        {
            if (A[i][cur_j] != 0.0)
            {
                t = -A[i][cur_j];
                A[i][cur_j] = 0.0;
                
                for (int j=cur_j+1;j<n;j++) A[i][j] += t * A[cur_j][j];
                for (int j=0;j<nb;j++) B[i][j] += t * B[cur_j][j];
            }
        }
    }
    
    for (int cur_j=n-1;cur_j>=1;cur_j--)
    {
        for (int i=0;i<cur_j;i++)
        {
            for (int j=0;j<nb;j++)
            {
                B[i][j] -= A[i][cur_j] * B[cur_j][j];
            }
        }
    }
}

double** mat_exp(double **M, int n)
{
    int q = 6;
    
    double norm = 0.0;
    for (int i=0;i<n;i++)
    {
        double cur_sum = 0.0;
        for (int j=0;j<n;j++) cur_sum += fabs(M[i][j]);
        if (cur_sum > norm) norm = cur_sum;
    }
    
    int s = max(0, int(log2(norm)) + 1);
    double t = pow(2.0, -s);
    
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            M[i][j] *= t;
        }
    }
    
    double **tmp = new double*[n];
    for (int i=0;i<n;i++) tmp[i] = new double[n];
    
    double **x = new double*[n];
    for (int i=0;i<n;i++)
    {
        x[i] = new double[n];
        for (int j=0;j<n;j++)
        {
            x[i][j] = M[i][j];
        }
    }
    
    double c = 0.5;
    
    double **d = new double*[n];
    double **ret = new double*[n];
    for (int i=0;i<n;i++)
    {
        ret[i] = new double[n];
        d[i] = new double[n];
        
        for (int j=0;j<n;j++)
        {
            ret[i][j] = (i == j) + c * M[i][j];
            d[i][j] = (i == j) - c * M[i][j];
        }
    }
    
    bool p = true;
    
    for (int k=2;k<=q;k++)
    {
        c *= (1.0 * (q - k + 1)) / (1.0 * (k * ((q << 1) - k + 1)));
        
        for (int i=0;i<n;i++)
        {
            for (int j=0;j<n;j++)
            {
                tmp[i][j] = 0.0;
                for (int k=0;k<n;k++)
                {
                    tmp[i][j] += M[i][k] * x[k][j];
                }
            }
        }
        
        for (int i=0;i<n;i++)
        {
            for (int j=0;j<n;j++)
            {
                x[i][j] = tmp[i][j];
                ret[i][j] += c * x[i][j];
                if (p) d[i][j] += c * x[i][j];
                else d[i][j] -= c * x[i][j];
            }
        }
        
        p = !p;
    }
    
    solve(d, ret, n, n);
    
    while (s--)
    {
        for (int i=0;i<n;i++)
        {
            for (int j=0;j<n;j++)
            {
                tmp[i][j] = 0.0;
                for (int k=0;k<n;k++)
                {
                    tmp[i][j] += ret[i][k] * ret[k][j];
                }
            }
        }
        
        for (int i=0;i<n;i++)
        {
            for (int j=0;j<n;j++)
            {
                ret[i][j] = tmp[i][j];
            }
        }
    }
    
    for (int i=0;i<n;i++)
    {
        delete[] tmp[i];
        delete[] x[i];
        delete[] d[i];
    }
    
    delete[] tmp;
    delete[] x;
    delete[] d;
    
    return ret;
}

int main()
{
    int len;
    double **mat;
    
    printf("Input matrix size: ");
    scanf("%d", &len);
    
    mat = new double*[len];
    
    printf("Input matrix:\n");
    for (int i=0;i<len;i++)
    {
        mat[i] = new double[len];
        for (int j=0;j<len;j++)
        {
            scanf("%lf", &mat[i][j]);
        }
    }
    
    printf("EXP(M) = \n");
    double **fin = mat_exp(mat, len);
    for (int i=0;i<len;i++)
    {
        for (int j=0;j<len;j++)
        {
            printf("%lf ", fin[i][j]);
        }
        printf("\n");
    }
    
    return 0;
}