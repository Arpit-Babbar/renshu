#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include <cassert>
using namespace std;

template<typename Pde_Solver>
void run_and_get_output(double n_points, double max_refinements,
                        Pde_Solver &solver); //Never forget to put & here!!