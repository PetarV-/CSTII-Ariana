#ifndef NSGAII_OBJ
#define NSGAII_OBJ

#include <vector>

typedef double func_t (std::vector<double>);
typedef func_t* pfunc_t;

std::vector<pfunc_t> get_objectives();

#endif