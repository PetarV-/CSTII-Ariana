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

#include <nsga2.h>

#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;
typedef long long lld;
typedef unsigned long long llu;

extern vector<function<double(vector<double>)> > get_objectives();

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        printf("Usage: ./nsga2_tester <input_parameter_file> <results_file>\n");
        return -1;
    }
    
    printf("Preparing input parameters...\n");
    
    // Prepare the input parameters
    FILE *f = fopen(argv[1], "r");
    
    nsga2_params params;
    fscanf(f, "%d", &params.pop_size);
    fscanf(f, "%d", &params.ft_size);
    fscanf(f, "%d", &params.obj_size);
    fscanf(f, "%d", &params.generations);
    fscanf(f, "%lf", &params.p_crossover);
    fscanf(f, "%lf", &params.p_mutation);
    fscanf(f, "%lf", &params.di_crossover);
    fscanf(f, "%lf", &params.di_mutation);
    
    params.var_lims.resize(params.ft_size);
    for (int i=0;i<params.ft_size;i++)
    {
        double x, y;
        fscanf(f, "%lf%lf", &x, &y);
        params.var_lims[i] = make_pair(x, y);
    }
    
    fclose(f);
    
    
    printf("Running NSGA-II...\n");
    
    // Run the algorithm
    NSGAII nsga2;
    vector<function<double(vector<double>)> > objs = get_objectives();
    vector<chromosome> final_generation = nsga2.optimise(params, objs);
    
    printf("Extracting Pareto front...\n");
    
    FILE *g = fopen(argv[2], "w");
    
    fprintf(g, "Pareto front:\n");
    
    int ii = 0;
    
    list<chromosome> pareto_front = nsga2.find_nondominated_front(final_generation);
    
    for (list<chromosome>::iterator it = pareto_front.begin(); it != pareto_front.end(); it++)
    {
        chromosome cur = (*it);
        fprintf(g, "Chromosome id %d: (", ii++);
        for (int i=0;i<nsga2.get_ft_size();i++) fprintf(g, (i == nsga2.get_ft_size() - 1) ? "%lf) -> (" : "%lf, ", cur.features[i]);
        for (int i=0;i<nsga2.get_obj_size();i++) fprintf(g, (i == nsga2.get_obj_size() - 1) ? "%lf)\n" : "%lf, ", cur.values[i]);
    }
    
    fprintf(g, "===========================\n");
    
    for (uint jj=ii;jj<final_generation.size();jj++)
    {
        chromosome cur = final_generation[jj];
        fprintf(g, "Chromosome id %d: (", jj++);
        for (int i=0;i<nsga2.get_ft_size();i++) fprintf(g, (i == nsga2.get_ft_size() - 1) ? "%lf) -> (" : "%lf, ", cur.features[i]);
        for (int i=0;i<nsga2.get_obj_size();i++) fprintf(g, (i == nsga2.get_obj_size() - 1) ? "%lf)\n" : "%lf, ", cur.values[i]);
    }
    
    printf("Pareto front extracted. Results saved in %s.\n", argv[2]);
    
    fclose(g);
    
    return 0;
}
