#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include "finite_volume_solver.h"
#include <cassert>

using namespace std;
Finite_Volume_Solver_1d::Finite_Volume_Solver_1d(double n_points, double cfl, string scheme,
                                           double running_time,
                                           string initial_data_indicator):
                                           n_points(n_points), cfl(cfl),
                                           running_time(running_time), scheme(scheme),
                                           initial_data_indicator(initial_data_indicator)
{
  coefficient = 1.0;
  sigma = cfl * coefficient / abs(coefficient); //sigma = coefficient*dt/h
  x_min = -1.0, x_max = 1.0;
  initial_data_function.set_initial_data(x_min,x_max,initial_data_indicator);
  //This would set our initial data function, inputted by user. 
  //The function changes with the domain, so we need x_min,x_max.
  h = (x_max - x_min) / (n_points);
  t = 0.0;
  dt = cfl * h / abs(coefficient);
  cout << "h = " << h << endl;
  cout << "dt = " << dt << endl;
  grid.resize(n_points);
  error.resize(n_points);
  solution_old.resize(n_points);
  solution.resize(n_points);
  solution_exact.resize(n_points);
  temp.resize(n_points);
}

void Finite_Volume_Solver_1d::make_grid()
{
    //For n = n_points, Grid = [x_0,x_1,..., x_{n-1}]
    //x_0 = x_min + 0.5*h
    //x_{j+1} - x_j = h = (x_max - x_min)/n
    //x_{n-1} = x_max - 0.5*h
    for (unsigned int i = 0; i < n_points; i++)
    {
        grid[i] = x_min + 0.5*h + i * h;
    }
}
void Finite_Volume_Solver_1d::refine()
{
  n_points = 2.0*n_points;
  sigma = cfl * coefficient / abs(coefficient); //sigma = coefficient*dt/h
  //This would set our initial data function, inputted by user. 
  //The function changes with the domain, so we need x_min,x_max.
  h = (x_max - x_min) / (n_points);
  t = 0.0;
  dt = cfl * h / abs(coefficient);
  cout << "Values after refinement are"<<endl;
  cout << "h = " << h << endl;
  cout << "dt = " << dt << endl;
  grid.resize(n_points);
  error.resize(n_points);
  solution_old.resize(n_points);
  solution.resize(n_points);
  solution_exact.resize(n_points);
  temp.resize(n_points);
}

void Finite_Volume_Solver_1d::set_initial_solution()
{
    for (unsigned int i = 0; i < n_points; i++)
    {
        solution[i] = initial_data_function.value(grid[i]);
    }
}

void Finite_Volume_Solver_1d::lax_wendroff()
{
    solution[0] = solution_old[0] 
                     - 0.5 * cfl * (solution_old[1] - solution_old[n_points - 1])
                     + 0.5 * cfl * cfl * (solution_old[n_points - 1] 
                                          - 2.0 * solution_old[0] + solution_old[1]);
    for (unsigned int j = 1; j < n_points - 1; ++j) //Loop over grid points
    {
        solution[j] = solution_old[j] 
                          - 0.5 * cfl * (solution_old[j + 1] - solution_old[j - 1]) 
                          + 0.5 * cfl * cfl * (solution_old[j - 1] 
                                               - 2.0 * solution_old[j] + solution_old[j + 1]);
    }
    solution[n_points - 1] = solution_old[n_points - 1] 
                                 - 0.5 * cfl * (solution_old[0] - solution_old[n_points - 2])
                                 + 0.5 * cfl * cfl * (solution_old[n_points - 2] 
                                                     - 2 * solution_old[n_points - 1]
                                                     + solution_old[0]);
}


void Finite_Volume_Solver_1d::evaluate_error_and_output_solution(int time_step_number)
{
    for (unsigned int i = 0; i < n_points; i++)
    {
        solution_exact[i] = initial_data_function.value(grid[i] - coefficient * t);
    }
    
    for (unsigned int j = 0; j < n_points; j++)
    {
        error[j] = abs(solution[j] - solution_exact[j]);
    }
    string solution_file_name = "solution_";
    solution_file_name += to_string(time_step_number) + ".txt";
    output_vectors_to_file(solution_file_name, grid, solution, solution_exact);
}

void Finite_Volume_Solver_1d::output_final_error()
{
    string error_file_name = "finalerror_";
    error_file_name += scheme + ".txt";
    output_vectors_to_file(error_file_name, grid, error);
}

void Finite_Volume_Solver_1d::get_error(vector<double> &l1_vector, 
                                     vector<double> &l2_vector,
                                     vector<double> &linfty_vector)
{
    double l1 = 0.0;
    double l2 = 0.0;
    double linfty = 0.0;
    for (int j = 0; j < n_points; j++) 
    {
        l1 = l1 + error[j]  * h;           // L1 error
        l2 = l2 + error[j] * error[j]  * h;// L2 error
        linfty = max(linfty, error[j]);            // L_infty error
    }
    l2 = sqrt(l2);
    l1_vector.push_back(l1);
    l2_vector.push_back(l2);
    linfty_vector.push_back(linfty);
}

