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

#include "multiplex.h"
#include "../matrix_lib/matrix_lib.h"

#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;
typedef long long lld;
typedef unsigned long long llu;

Multiplex::Multiplex(int n, int L) : n(n), L(L)
{
    this -> layers.resize(L);
    for (int i=0;i<L;i++)
    {
        layers[i] = new SimpleGraphLayer(n);
    }
    
    M = new double***[L];
    for (int i=0;i<L;i++)
    {
        M[i] = new double**[L];
        for (int j=0;j<L;j++)
        {
            M[i][j] = new double*[n];
            for (int k=0;k<n;k++)
            {
                M[i][j][k] = new double[n];
                for (int l=0;l<n;l++)
                {
                    M[i][j][k][l] = 0.0;
                }
            }
        }
    }
}

Multiplex::Multiplex(int n, int L, double ****A) : n(n), L(L)
{
    this -> layers.resize(L);
    for (int i=0;i<L;i++)
    {
        layers[i] = new SimpleGraphLayer(n, A[i][i]);
    }
    
    M = new double***[L];
    for (int i=0;i<L;i++)
    {
        M[i] = new double**[L];
        for (int j=0;j<L;j++)
        {
            M[i][j] = new double*[n];
            for (int k=0;k<n;k++)
            {
                M[i][j][k] = new double[n];
                for (int l=0;l<n;l++)
                {
                    M[i][j][k][l] = A[i][j][k][l];
                }
            }
        }
    }
}

Multiplex::Multiplex(int n, vector<AbstractGraphLayer*> lyrs) : n(n), L(lyrs.size())
{
    this -> layers.resize(L);
    
    for (int i=0;i<L;i++)
    {
        assert(lyrs[i] -> get_n() == n);
        this -> layers[i] = lyrs[i];
    }
    
    M = new double***[L];
    for (int i=0;i<L;i++)
    {
        M[i] = new double**[L];
        for (int j=0;j<L;j++)
        {
            M[i][j] = new double*[n];
            for (int k=0;k<n;k++)
            {
                M[i][j][k] = new double[n];
                for (int l=0;l<n;l++)
                {
                    if (i != j) M[i][j][k][l] = 0.0;
                    else M[i][j][k][l] = lyrs[i] -> get_adj(k, l);
                }
            }
        }
    }
}


Multiplex::Multiplex(int n, std::vector<AbstractGraphLayer*> lyrs, double **omega) : n(n), L(lyrs.size())
{
    this -> layers.resize(L);
    
    for (int i=0;i<L;i++)
    {
        assert(lyrs[i] -> get_n() == n);
        this -> layers[i] = lyrs[i];
    }
    
    M = new double***[L];
    for (int i=0;i<L;i++)
    {
        M[i] = new double**[L];
        for (int j=0;j<L;j++)
        {
            M[i][j] = new double*[n];
            for (int k=0;k<n;k++)
            {
                M[i][j][k] = new double[n];
                for (int l=0;l<n;l++)
                {
                    if (i != j)
                    {
                        if (k != l) M[i][j][k][l] = 0.0;
                        else M[i][j][k][l] = omega[i][j];
                    }
                    else M[i][j][k][l] = lyrs[i] -> get_adj(k, l);
                }
            }
        }
    }
}

Multiplex::Multiplex(int n, std::vector<AbstractGraphLayer*> lyrs, double ****inter_lyr) : n(n), L(lyrs.size())
{
    this -> layers.resize(L);
    
    for (int i=0;i<L;i++)
    {
        assert(lyrs[i] -> get_n() == n);
        this -> layers[i] = lyrs[i];
    }
    
    M = new double***[L];
    for (int i=0;i<L;i++)
    {
        M[i] = new double**[L];
        for (int j=0;j<L;j++)
        {
            M[i][j] = new double*[n];
            for (int k=0;k<n;k++)
            {
                M[i][j][k] = new double[n];
                for (int l=0;l<n;l++)
                {
                    if (i != j) M[i][j][k][l] = inter_lyr[i][j][k][l];
                    else M[i][j][k][l] = lyrs[i] -> get_adj(k, l);
                }
            }
        }
    }
}

double Multiplex::get_edge(int lyr1, int lyr2, int n1, int n2)
{
    return M[lyr1][lyr2][n1][n2];
}

void Multiplex::set_edge(int lyr1, int lyr2, int n1, int n2, double val)
{
    M[lyr1][lyr2][n1][n2] = val;
    if (lyr1 == lyr2) layers[lyr1] -> set_adj(n1, n2, val);
}

double** Multiplex::get_matrix_form()
{
    double** ret = new double*[n*L];
    for (int i=0;i<n*L;i++)
    {
        ret[i] = new double[n*L];
        
        int a = i / n;
        int c = i % n;
        
        for (int j=0;j<n*L;j++)
        {
            int b = j / n;
            int d = j % n;
            ret[i][j] = M[a][b][c][d];
        }
    }
    return ret;
}

double**** Multiplex::get_communicability_matrix()
{
    double **mat = get_matrix_form();
    double **ex = mat_exp(mat, n*L);
    
    double ****ret = new double***[L];
    for (int i=0;i<L;i++)
    {
        ret[i] = new double**[L];
        for (int j=0;j<L;j++)
        {
            ret[i][j] = new double*[n];
            for (int k=0;k<n;k++)
            {
                ret[i][j][k] = new double[n];
                for (int l=0;l<n;l++)
                {
                    ret[i][j][k][l] = mat[i*L+k][j*L+l];
                }
            }
        }
    }
    
    for (int i=0;i<n*L;i++)
    {
        delete[] mat[i];
        delete[] ex[i];
    }
    
    delete[] mat;
    delete[] ex;
    
    return ret;
}