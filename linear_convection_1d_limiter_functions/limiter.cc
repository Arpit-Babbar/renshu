#include <cmath>
#include <iostream>
#include <fstream>

#include <cassert>
#include <string>
#include <stdio.h>
#include <sys/time.h>
#include "run_and_get_output.h"
#include "finite_volume_solver.h"
using namespace std;

//This function will output the flux f_{i+1/2} from input v_i,v_{i-1}

double minmod(double numerator,double denominator)
{
  if (numerator * denominator <= 0.0) 
      return 0.0;
  else 
      return min(numerator / denominator,1.0) * denominator;
}

double superbee(double numerator,double denominator)
{
  if (numerator * denominator <= 0.0) 
      return 0.0;
  else 
      {
      double r = numerator/denominator;
      double temp = max(0.0,min(2.0*r,1.0));
      return max(temp, min(r,2.0)) * denominator;
      }
}

//back_diff will mean uj - ujm1 for some j. Similar for cent_diff,fwd_diff.
double minmod3(double back_diff,double cent_diff,double fwd_diff)
{
  if (  (back_diff*cent_diff <= 0.0) ||
        (cent_diff*fwd_diff <= 0.0) )
      return 0.0;
  else //s min(|a|,|b|,|c|) where s = sign(a)=sign(b)=sign(c)
      {
        double temp = min(abs(back_diff),abs(cent_diff));
        temp = min(temp,abs(fwd_diff));
        return (back_diff/abs(back_diff)) * temp;
      }
}


class Limiter_1d : public Finite_Volume_Solver_1d
{
public:
    Limiter_1d(const double n_points, const double cfl,
               string method, const double running_time,
               string initial_data_indicator, string limiter);
    void run();
protected:
    void make_grid();
    void rhs_function();

    void foup(); //First order upwind scheme  
    void ssp_rk2_solver();
    void ssp_rk3_solver();

    double reconstructor(int i);

    void rhs_soup();

    string limiter;

    void soup_rk2(); //Second order upwind Scheme
    void soup_rk3();
};

Limiter_1d::Limiter_1d(const double n_points, const double cfl,
            string method, const double running_time,
            string initial_data_indicator, string limiter)
            :Finite_Volume_Solver_1d(n_points,cfl,
            method, running_time,initial_data_indicator), limiter(limiter){}

void Limiter_1d::make_grid()
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

void Limiter_1d::foup() //first order upwind scheme
{
  solution[0] = (1.0 - sigma) * solution_old[0] + cfl * solution_old[n_points - 1];
  for (int i = 1; i < n_points; i++)
  {
    solution[i] = (1.0 - sigma) * solution_old[i] + cfl * solution_old[i - 1];
  }
}

double Limiter_1d::reconstructor(int j) //Gives v_{i-1/2}^L = v_i + 0.5*Phi(i)
{
  double ujm2,ujm1,uj/*,ujp1*/;//u_{j-2},u_{j-1},u_j,u_{j+1} isolated
  //for j=n-1,0,1,2 cases.
  //Though it is wasteful in most cases, it adds to readability.
  uj = solution[j];
  /*if (j==n_points-1)
    ujp1 = solution[0];
  else 
    ujp1 = solution[j+1];*/ //ujp1 is not needed in our limiters.
  if (j==0)
    {
      ujm2 = solution[n_points-2];
      ujm1 = solution[n_points-1];
    }
  else if (j==1)
    {
      ujm2 = solution[n_points-1];
      ujm1 = solution[j-1];
    }
  else 
    {
      ujm2 = solution[j-2];
      ujm1 = solution[j-1];
    }
    
  if (limiter == "none")
  {
    return ujm1 + 0.5*(ujm1-ujm2);
  }
  else if (limiter == "minmod")
  {
    return ujm1 + 0.5*minmod(uj - ujm1,ujm1 -ujm2);
  }
  else if (limiter == "superbee")
  {
    return ujm1 + 0.5*superbee(uj-ujm1,ujm1-ujm2);
  }
  else if (limiter == "vanleer")
  {
    double beta = 2.0;
    return ujm1 + 0.5 * minmod3(beta*(ujm1-ujm2),
                                0.5*(uj-ujm2),beta*(uj-ujm1));
  }
  else
  {
    cout << "Incorrect limiter inputted"<<endl;
    assert(false);
  }
}

void Limiter_1d::rhs_function()
{
  fill((*rhs).begin(),(*rhs).end(),0.0); //Sets (*rhs) vector to zero.
  //There are n_points+1 points on which the flux needs to be evaluated
  //But, the last and first flux are the same, so we'd not evaluate them
  
  //Also, each flux finds its use in computing the RHS at two of the 
  //'FVM grid' points, so we'd put it in both the places in same iteration.
  double flux;
  //Recall that we had defined flux_{j+1/2} = coefficient*u_{j+1/2}^L.
  //We'd just compute that u_{j+1/2}^L, and store it as ul.
  double ul;
  //For usual reasons, we need to consider the first and last flux separately
  //First we compute f_{-1/2}
  ul = reconstructor(0);//u_{0-1/2}^L = u_{-1/2}^L
  flux = coefficient*ul; //f_{-1/2}
  (*rhs)[0]          +=  flux / h; //rhs[0] = -(f_{1/2} - f_{-1/2})/h,
  (*rhs)[n_points-1] += -flux / h;
  //rhs[n-1]=-(f_{n-1/2}-f_{n-3/2})=-(f_{-1/2}-f_{n-3/2})
  
  //Next compute f_{1/2}
  ul = reconstructor(1); //u_{1-1/2}^L=u_{1/2}^L
  flux = coefficient*ul; //f_{1/2}
  (*rhs)[1] +=  flux/h; //rhs[1] =-(f_{3/2}-f_{1/2})/h
  (*rhs)[0] += -flux/h; //rhs[0] =-(f_{1/2}-f_{-1/2})/h
  for (int i = 2; i<n_points;i++) //Only n_points-1 fluxes to be computed
  {
    ul = reconstructor(i); //u_{i-1/2}^L
    flux = coefficient * ul;//f_{i-1/2}
    (*rhs)[i]  +=  flux/h; //rhs[i]  =-(f_{i+1/2}-f_{i-1/2})/h
    (*rhs)[i-1]+= -flux/h; //rhs[i-1]=-(f_{i-1/2}-f_{i-3/2})/h
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
      solution_old = solution;//update solution_old to be used at next time step
      if (method == "lw")
        lax_wendroff();
      else if (method == "foup")
        foup();
      else if (method == "soup_rk3")
        ssp_rk3_solver();
      else if (method == "soup_rk2")
        ssp_rk2_solver();
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
      cout << "Incorrect format, use" << endl;
      cout << "./output method cfl running_time initial_data n_refinements limiter"<<endl;
      cout << "Choice for method - lw,foup,soup_rk3,soup_rk2 ." << endl;
      cout << "Choices for initial_data - smooth_sine,hat,step,cts_sine . " <<endl;
      cout << "Optional choices for limiter are - minmod,superbee"<<endl;
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
  string limiter;
  if (argc == 6)
    limiter = "none";
  else
    limiter = argv[6];
  cout << "Chosen limiter is "<< limiter;
  Limiter_1d solver(n_points, cfl, method, running_time,
                                  initial_data_indicator,limiter);
  run_and_get_output(n_points, max_refinements, solver);
}
