#include <cmath>
#include <iostream>
#include <fstream>

#include <cassert>
#include <string>
#include <stdio.h>
#include <sys/time.h>
#include "vector_add.h"
#include "run_and_get_output.h"
#include "finite_volume_solver.h"
using namespace std;

//This function will output the flux f_{i+1/2} from input v_i,v_{i-1}

class Limiter_1d : public Finite_Volume_Solver_1d
{
public:
    Limiter_1d(const double n_points, const double cfl,
               string method, const double running_time,
               string initial_data_indicator)
               :Finite_Volume_Solver_1d(n_points,cfl,
               method, running_time,initial_data_indicator){}
    void run();
protected:
    void rhs_function();

    void foup(); //First order upwind scheme  
    void ssp_rk2_solver();
    void ssp_rk3_solver();
};


void Limiter_1d::foup() //first order upwind scheme
{
  solution[0] = (1.0 - sigma) * solution_old[0] + cfl * solution_old[n_points - 1];
  for (int i = 1; i < n_points; i++)
  {
    solution[i] = (1.0 - sigma) * solution_old[i] + cfl * solution_old[i - 1];
  }
}

void Limiter_1d::rhs_function()
{
  fill((*rhs).begin(),(*rhs).end(),0.0); //Sets (*rhs) vector to zero.
  //There are n_points+1 points on which the flux needs to be evaluated
  //But, the last and first flux are the same, so we'd not evaluate them
  
  //Also, each flux, other than first and last, 
  //finds its use in computing the RHS at two of the 
  //'FVM grid' points, so we'd put it in both the places in same iteration.
  
  //The first and last flux points are required by first and last grid 
  //points, so only two points require them in total.

  double flux;
  //Recall that we had defined flux_{j+1/2} = coefficient*u_{j+1/2}^L.
  //We'd just compute that u_{j+1/2}^L, and store it as ul.
  double ul;
  //For usual reasons, we need to consider the first two fluxes separately
  //First we compute f_{-1/2}
  ul = solution[n_points-1] + 0.5 *(solution[n_points-1]-solution[n_points-2]);
  flux = coefficient*ul; //f_{-1/2}
  (*rhs)[0]          +=  (flux/h); //rhs[0] = -(f_{1/2} - f_{-1/2})/h,
  (*rhs)[n_points-1] -=  (flux/h);
  //rhs[n-1]=-(f_{n-1/2}-f_{n-3/2})=-(f_{-1/2}-f_{n-3/2})
  
  //Next compute f_{1/2}
  ul = solution[0] + 0.5*(solution[0] - solution[n_points-1]); //u_{1/2}^L
  flux = coefficient*ul; //f_{1/2}
  (*rhs)[1] +=  (flux/h); //rhs[1] =-(f_{3/2}-f_{1/2})/h
  (*rhs)[0] -=  (flux/h); //rhs[0] =-(f_{1/2}-f_{-1/2})/h
  
  //Now computing f_{i-1/2}
  for (int i = 2; i<n_points;i++) //Only n_points-1 fluxes to be computed
  {
    ul = solution[i-1] + 0.5 *(solution[i-1]-solution[i-2]); //u_{i-1/2}^L
    flux = coefficient * ul;//f_{i-1/2}
    (*rhs)[i]  +=  (flux/h); //rhs[i]  =-(f_{i+1/2}-f_{i-1/2})/h
    (*rhs)[i-1]-=  (flux/h); //rhs[i-1]=-(f_{i-1/2}-f_{i-3/2})/h
  }
}

void Limiter_1d::ssp_rk3_solver()
{
    rhs = &temp;
    rhs_function(); //Put temp = rhs(solution)=rhs(solution_old)
    add(solution_old,dt,temp,solution); //Put solution = solution_old + dt*rhs(solution_old)
    rhs_function(); //Put temp = f(y)
              //We want to put k2 = 3/4 * U^n + 1/4*(k1 + dt*rhs(k1))
    add(solution,dt,temp,temp);
    add(0.75,solution_old,0.25,temp, solution);
    rhs_function(); //r = f(y)
    add(solution,dt,temp,solution);
    add(1.0/3.0, solution_old, 2.0/3.0, solution, solution);
}

void Limiter_1d::ssp_rk2_solver()
{
    rhs = &temp;
    rhs_function();
    add(solution_old,dt,temp,solution);
    rhs_function();
    add(solution,dt,temp,temp);
    add(0.5,solution_old,0.5,temp,solution);    
}

void Limiter_1d::run()
{
    make_grid();
    set_initial_solution(); //sets solution to be the initial data
    int time_step_number = 0;
    evaluate_error_and_output_solution(time_step_number);
    while (t < running_time) //compute solution at next time step using solution_old
    {
        solution_old = solution;//update solution_old to be used at this time step
        if (method == "lw")
            lax_wendroff();
        else if (method == "foup")
            foup();
        else if (method == "soup_rk3" ||  method == "soup_rk3_minmod" || 
                 method == "soup_rk3_superbee")
            ssp_rk3_solver();
            //rk3_solver in which the rhs would depend on method
        else if (method == "soup_rk2" || method == "soup_rk2_minmod" ||
                 method == "soup_rk2_superbee")
            ssp_rk2_solver();
            //rk2_solver in which the rhs would depend on method
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
    if (argc != 6)
    {
        cout << "Incorrect format, use" << endl;
        cout << "./output lw/rk4/rk3/rk3_new/rk2(for the respective method) cfl running_time 0/1/2(for " ;
        cout << "sin/hat/discts initial data) max_refinements" << endl;
        assert(false);
    }
    string method = argv[1];
    cout << "method = " << method << endl;
    double n_points = 50.0; //Number of x_j type points, not x_{j+1/2} type.
    double cfl = stod(argv[2]);
    cout << "cfl = " << cfl << endl;
    double running_time = stod(argv[3]);
    cout << "running_time = " << running_time << endl;
    string initial_data_indicator = argv[4];
    cout << "initial_data_indicator = " << initial_data_indicator << endl;
    unsigned int max_refinements = stoi(argv[5]);
    cout << "max_refinements = " << max_refinements <<endl;
    Limiter_1d solver(n_points, cfl, method, running_time,
                                    initial_data_indicator);
    run_and_get_output(n_points, cfl, method, running_time, initial_data_indicator,
                       max_refinements, solver);
}
