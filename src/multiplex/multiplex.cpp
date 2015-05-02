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

#include <multiplex.h>
#include <matrix_lib.h>
#include <nsga2.h>
#include <vector_cmp.h>

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

Multiplex::Multiplex(int n, vector<AbstractGraphLayer*> &lyrs) : n(n), L(lyrs.size())
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


Multiplex::Multiplex(int n, vector<AbstractGraphLayer*> &lyrs, double **omega) : n(n), L(lyrs.size())
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

Multiplex::Multiplex(int n, vector<AbstractGraphLayer*> &lyrs, double ****inter_lyr) : n(n), L(lyrs.size())
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

Multiplex::~Multiplex()
{
    for (int i=0;i<L;i++) delete layers[i];
    for (int i=0;i<L;i++)
    {
        for (int j=0;j<L;j++)
        {
            for (int k=0;k<n;k++) delete[] M[i][j][k];
            delete[] M[i][j];
        }
        delete[] M[i];
    }
    delete[] M;
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

double** Multiplex::get_aggregate_matrix(bool normalise)
{
    double ****G = get_communicability_matrix();
    
    double **ret = new double*[n];
    
    for (int i=0;i<n;i++)
    {
        ret[i] = new double[n];
        for (int j=0;j<n;j++)
        {
            if (i == j) ret[i][j] = 0.0;
            else
            {
                bool cont_zero = false;
                double h_mean_den = 0.0;
                for (int k=0;k<L;k++)
                {
                    for (int l=0;l<L;l++)
                    {
                        h_mean_den += 1.0 / G[k][l][i][j];
                        if (G[k][l][i][j] == 0) cont_zero = true;
                    }
                }
                ret[i][j] = (cont_zero ? 0 : ((1.0 * L * L) / h_mean_den));
            }
        }
    }
    
    if (normalise)
    {
        for (int i=0;i<n;i++)
        {
            double sum = 0.0;
            for (int j=0;j<n;j++)
            {
                sum += ret[i][j];
            }
            if (sum >= 1e-6)
            {
                for (int j=0;j<n;j++)
                {
                    ret[i][j] /= sum;
                }
            }
        }
    }
    
    for (int i=0;i<L;i++)
    {
        for (int j=0;j<L;j++)
        {
            for (int k=0;k<n;k++)
            {
                delete[] G[i][j][k];
            }
            delete[] G[i][j];
        }
        delete[] G[i];
    }
    delete[] G;
    
    return ret;
}

void Multiplex::sync()
{
    for (int i=0;i<L;i++)
    {
        for (int j=0;j<n;j++)
        {
            for (int k=0;k<n;k++)
            {
                M[i][i][j][k] = layers[i] -> get_adj(j, k);
            }
        }
    }
}

void Multiplex::set_omega(double **omega, bool normalise, bool upsync)
{
    for (int i=0;i<L;i++)
    {
        for (int j=0;j<L;j++)
        {
            for (int k=0;k<n;k++)
            {
                for (int l=0;l<n;l++)
                {
                    if (i != j)
                    {
                        if (k != l) M[i][j][k][l] = 0.0;
                        else M[i][j][k][l] = omega[i][j];
                    }
                }
            }
        }
    }
    
    // now normalise if desired
    if (normalise)
    {
        for (int i=0;i<L;i++)
        {
            for (int k=0;k<n;k++)
            {
                double sum = 0.0;
                for (int j=0;j<L;j++)
                {
                    for (int l=0;l<n;l++)
                    {
                        sum += M[i][j][k][l];
                    }
                }
                for (int j=0;j<L;j++)
                {
                    for (int l=0;l<n;l++)
                    {
                        M[i][j][k][l] /= sum;
                        if (i == j && upsync) layers[i] -> set_adj(k, l, M[i][j][k][l]);
                    }
                }
            }
        }
    }
}

void Multiplex::train(vector<vector<vector<double> > > &train_set)
{
    // Train all the layers individually (as before)
    for (int l=0;l<L;l++)
    {
        vector<vector<double> > curr_set(train_set.size(), vector<double>(n));
        for (uint i=0;i<train_set.size();i++)
        {
            for (int j=0;j<n;j++)
            {
                curr_set[i][j] = train_set[i][j][l];
            }
        }
        layers[l] -> train(curr_set);
    }
    
    sync();
    
    // Define the lambdas that calculate likelihoods for a given omega
    objectives.resize(train_set.size());
    for (uint t=0;t<train_set.size();t++)
    {
        objectives[t] = [this, t, &train_set] (vector<double> X) -> double
        {
            int ii = 0;
            double **temp_omega = new double*[L];
            for (int i=0;i<L;i++)
            {
                temp_omega[i] = new double[L];
                for (int j=0;j<L;j++)
                {
                    if (i != j)
                    {
                        temp_omega[i][j] = X[ii++];
                    }
                    else temp_omega[i][j] = 0.0;
                }
            }
            
            set_omega(temp_omega, true, false);
            
            for (int i=0;i<L;i++) delete[] temp_omega[i];
            delete[] temp_omega;
            
            return -log_likelihood(train_set[t]);
        };
    }
    
    // Prepare the input parameters for NSGA-II
    string filename = "param.in";
    const int pop_size = 100;
    const int ft_size = L * (L - 1);
    const int obj_size = train_set.size();
    const int generations = 200;
    const double p_crossover = 0.9;
    const double p_mutation = 1.0;
    const double di_crossover = 10.0;
    const double di_mutation = 100.0;
    
    FILE *f = fopen(filename.c_str(), "w");
    fprintf(f, "%d\n", pop_size);
    fprintf(f, "%d\n", ft_size);
    fprintf(f, "%d\n", obj_size);
    fprintf(f, "%d\n", generations);
    fprintf(f, "%lf\n%lf\n", p_crossover, p_mutation);
    fprintf(f, "%lf\n%lf\n", di_crossover, di_mutation);
    for (int i=0;i<ft_size;i++)
    {
        fprintf(f, "0.000001 1.0\n");
    }
    
    fclose(f);
    
    // Run the algorithm
    NSGAII nsga2;
    vector<chromosome> candidates = nsga2.optimise((char*)filename.c_str(), objectives);
    
    
    // Evaluate the best choice of omega
    int best = -1;
    double min_sum = -1.0;
    for (uint i=0;i<candidates.size();i++)
    {
        sort(candidates[i].values.begin(), candidates[i].values.end());
        double curr_sum = 0.0;
        for (uint j=0;j<train_set.size();j++)
        {
            curr_sum += (j + 1) * candidates[i].values[j];
        }
        if (best == -1 || curr_sum < min_sum)
        {
            best = i;
            min_sum = curr_sum;
        }
    }
    
    // Adjust the parameters accordingly
    double **fin_omega = new double*[L];
    for (int i=0;i<L;i++)
    {
        fin_omega[i] = new double[L];
        for (int j=0;j<L;j++)
        {
            fin_omega[i][j] = candidates[best].features[i*L + j];
        }
    }
    
    set_omega(fin_omega, true, false);
    
    for (int i=0;i<L;i++) delete[] fin_omega[i];
    delete[] fin_omega;
    
    printf("TRANSITION MATRIX:\n");
    for (int i=0;i<L;i++)
    {
        for (int j=0;j<L;j++)
        {
            printf("%lf ", M[i][j][0][0]);
        }
        printf("\n");
    }
}

double Multiplex::log_likelihood(vector<vector<double> > &test_data)
{
    // First aggregate this multiplex
    double **aggr = this -> get_aggregate_matrix(true);
    
    // Then construct a representative multiplex of the test data
    vector<AbstractGraphLayer*> test_lyrs(L);
    for (int l=0;l<L;l++)
    {
        test_lyrs[l] = new SimpleGraphLayer(n);
        vector<vector<double> > curr_set(1, vector<double>(n));
        for (int j=0;j<n;j++)
        {
            curr_set[0][j] = test_data[j][l];
        }
        test_lyrs[l] -> train(curr_set);
    }
    
    double **omega = new double*[L];
    for (int i=0;i<L;i++)
    {
        omega[i] = new double[L];
        for (int j=0;j<L;j++)
        {
            omega[i][j] = M[i][j][0][0];
        }
    }
    
    Multiplex* test_mux = new Multiplex(n, test_lyrs, omega);
    
    // Aggregate it
    double **aggr_test = test_mux -> get_aggregate_matrix(true);
    
    // Then calculate the distance of those two matrices
    double ret = 0.0;
    
    for (int i=0;i<n;i++)
    {
        for (int j=0;j<n;j++)
        {
            ret += (aggr[i][j] - aggr_test[i][j]) * (aggr[i][j] - aggr_test[i][j]);
        }
    }
    
    for (int i=0;i<n;i++) delete[] aggr[i];
    delete[] aggr;
    
    for (int i=0;i<L;i++) delete[] omega[i];
    delete[] omega;
    
    delete test_mux;
    
    for (int i=0;i<n;i++) delete[] aggr_test[i];
    delete[] aggr_test;
    
    // As the desired result is similarity rather than distance,
    // return the negative log
    return -log(sqrt(ret));
}

void Multiplex::dump_muxviz_data(char *nodes_filename, char *base_layers_filename)
{
    FILE *f = fopen(nodes_filename, "w");
    
    fprintf(f, "nodeID\n");
    
    for (int i=1;i<=n;i++)
    {
        fprintf(f, "%d\n", i);
    }
    
    fclose(f);
    
    printf("Node data successfully written to %s.\n", nodes_filename);
    
    for (int i=0;i<L;i++)
    {
        char curr_lyr_filename[150];
        sprintf(curr_lyr_filename, "%s_%d", base_layers_filename, i+1);
        FILE *g = fopen(curr_lyr_filename, "w");
        for (int j=0;j<n;j++)
        {
            for (int k=0;k<n;k++)
            {
                fprintf(g, "%d %d %lf\n", j+1, k+1, layers[i] -> get_adj(j, k));
            }
        }
        fclose(g);
        printf("Layer %d data successfully written to %s.\n", i+1, curr_lyr_filename);
    }
    
    printf("Done.\n");
}
