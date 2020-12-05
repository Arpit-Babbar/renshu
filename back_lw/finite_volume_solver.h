#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "initial_conditions.h"
#include "vector_upgrade.h"
#include <cassert>

using namespace std;

class Finite_Volume_Solver_1d
{
public:
    Finite_Volume_Solver_1d(const double n_points, const double cfl,
                         string scheme, const double running_time,
                         string initial_data_indicator);
    //and which scheme to use - Lax-Wendroff or RK4
   //true when output is to be given, and false when it doesn't;
    void refine();
    void output_final_error();
    void get_error(vector<double> &l1_vector, vector<double> &l2_vector,
                   vector<double> &linfty_vector);
protected:
    void make_grid();
    void set_initial_solution();
    //removed initial_data vector, as it was wasteful
    //We'd use this to compute the rk4 slopes k_i's and temporarily store them in solution_new.
    //More precisely, it does solution = solution_old + factor * u
    void lax_wendroff();  
    
    void rhs_function();
    
    void evaluate_error_and_output_solution(const int time_step_number);

    double coefficient, x_min, x_max;

    vector<double> grid;

    vector<double> solution_old; //Solution at previous step
    vector<double> solution; //Solution at present step

    vector<double> *rhs; vector<double> temp;

    vector<double> solution_exact; //Exact solution at present time step

    vector<double> error; //This will store the maximum error at a grid point in all time-steps
    double l2_error;

    //h denotes the spatial distance between grid points
    double n_points, h, dt,t, cfl, sigma, running_time; //h = 1/n_points just included for easy typing
    //sigma = a h/dt = cfl * sign(a)

    string scheme;
    string initial_data_indicator;
    Initial_Data initial_data_function;
};
