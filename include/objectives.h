#ifndef NSGAII_OBJ
#define NSGAII_OBJ

#include <functional>
#include <vector>

std::vector<std::function<double(std::vector<double>)> > get_objectives();

#endif
