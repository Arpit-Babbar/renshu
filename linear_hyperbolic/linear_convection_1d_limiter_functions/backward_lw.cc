#include <cmath>
#include <iostream>
#include <fstream>

#include <cassert>
#include <string>
#include <stdio.h>
#include <sys/time.h>
#include "finite_volume_solver.h"
#include "run_and_get_output.h"
using namespace std;


class Solver : public Finite_Volume_Solver_1d
{
public:
    //and which scheme to use - Lax-Wendroff or RK4
    Solver(const double n_points, const double cfl,
               string scheme, const double running_time,
               string initial_data_indicator);
    void run();
protected:
    void make_grid();
    void rhs_function();
    void back_lw();
};

Solver::Solver(const double n_points, const double cfl,
               string scheme, const double running_time,
               string initial_data_indicator)
               :Finite_Volume_Solver_1d(n_points,cfl,
               scheme, running_time,initial_data_indicator){}

void Solver::make_grid()
{
    //Finite difference grid
    for (unsigned int i = 0; i < n_points; i++)
    {
        grid[i] = x_min + i * h; //grid that doesn't include x_max.
    }
}

void Solver::back_lw()
{
  solution[0] = (1.0 - 1.5 * sigma + 0.5 * sigma * sigma) * solution_old[0] 
                  + (2.0 * sigma - sigma * sigma) * solution_old[n_points-1]
                  + (-0.5 * sigma - 0.5 * sigma * sigma) * solution_old[n_points-2] ;
  solution[1] = (1.0 - 1.5 * sigma + 0.5 * sigma * sigma) * solution_old[1]
                  + (2.0 * sigma - sigma * sigma) * solution_old[0]
                  + (-0.5 * sigma - 0.5 * sigma * sigma) * solution_old[n_points-1];
  for (int i = 2; i<n_points; i++)
  {
    solution[i] = (1.0 - 1.5 * sigma + 0.5 * sigma * sigma) * solution_old[i]
                  + (2.0 * sigma - sigma * sigma) * solution_old[i-1]
                  + (-0.5 * sigma - 0.5 * sigma * sigma) * solution_old[i-2];
  }
}

void Solver::run()
{
    make_grid();
    set_initial_solution(); //sets solution to be the initial data
    int time_step_number = 0;
    evaluate_error_and_output_solution(time_step_number);
    while (t < running_time) //compute solution at next time step using solution_old
    {
      solution_old = solution;//update solution_old to be used at next time step
      if (scheme == "back_lw")
        back_lw();
      else
        {
          cout << "Incorrect scheme chosen "<<endl;
          assert(false);
        }
      time_step_number += 1;
      t = t + dt; 
      evaluate_error_and_output_solution(time_step_number);
    }
    cout << "For n_points = " << n_points<<", we took ";
    cout << time_step_number << " steps." << endl;
}


int main(int argc, char **argv)
{
  if ((argc != 6) && (argc != 7))
  {
      cout << "Incorrect number of arguments, use format" << endl;
      cout << "./output scheme cfl running_time initial_data n_refinements limiter"<<endl;
      cout << "Choices for scheme - lw,foup,soup_rk3,soup_rk2 ." << endl;
      cout << "Choices for initial_data - smooth_sine,hat,step,cts_sine . " <<endl;
      cout << "Blank limiter slot would run the scheme without a limiter "<<endl;
      cout << "Limiter choices are minmod,superbee,vanleer"<<endl;
      assert(false);
  }
  string scheme = argv[1];
  cout << "scheme = " << scheme << endl;
  double n_points = 50.0; //Number of x_j type points, not x_{j+1/2} type.
  double cfl = stod(argv[2]);
  cout << "cfl = " << cfl << endl;
  double running_time = stod(argv[3]);
  cout << "running_time = " << running_time << endl;
  string initial_data_indicator = argv[4];
  cout << "initial_data_indicator = " << initial_data_indicator << endl;
  unsigned int max_refinements = stoi(argv[5]);
  cout << "max_refinements = " << max_refinements <<endl;
  string limiter;
  if (argc == 6)
    limiter = "none";
  else
    limiter = argv[6];
  cout << "Chosen limiter is "<< limiter<< endl;
  Solver solver(n_points, cfl, scheme, running_time,
                                  initial_data_indicator);
  run_and_get_output(n_points, max_refinements, solver);
}

