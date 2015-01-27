/*
 Petar 'PetarV' Velickovic
 Algorithm: NSGA-II
*/

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
#include <chrono>
#include <random>

#include "objectives.h"

#define EPS 1e-9
#define INF 987654321

#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;
typedef unsigned int uint;
typedef long long lld;
typedef unsigned long long llu;

struct chromosome
{
    int rank;
    double distance;
    int sort_key; // for sorting purposes
    vector<double> features;
    vector<double> values;
    
    bool operator <(const chromosome &c) const
    {
        if (rank != c.rank) return (rank < c.rank);
        else return (distance > c.distance);
    }
};

list<chromosome> find_nondominated_front(vector<chromosome> &P);

vector<vector<chromosome> > fast_nondominated_sort(vector<chromosome> &P);

void crowding_distance_assignment(vector<chromosome> &I);

vector<chromosome> run(char *input_parameter_file);
