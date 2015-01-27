/*
 Petar 'PetarV' Velickovic
 Algorithm: NSGA-II
*/

#ifndef NSGAII
#define NSGAII

#include <vector>

struct chromosome
{
    int rank;
    double distance;
    int sort_key; // for sorting purposes
    std::vector<double> features;
    std::vector<double> values;
    
    bool operator <(const chromosome &c) const
    {
        if (rank != c.rank) return (rank < c.rank);
        else return (distance > c.distance);
    }
};

int get_ft_size();
int get_obj_size();

std::list<chromosome> find_nondominated_front(std::vector<chromosome> &P);

std::vector<std::vector<chromosome> > fast_nondominated_sort(std::vector<chromosome> &P);

void crowding_distance_assignment(std::vector<chromosome> &I);

std::vector<chromosome> optimise(char *input_parameter_file);

#endif
