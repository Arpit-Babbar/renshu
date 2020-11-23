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
double soup_reconstructor(double v1, double v0) //v_i, v_{i-1}
{
  return v1 + 0.5 * (v1 - v0);
}

class Limiter_1d : public Finite_Volume_Solver_1d
{
public:
    //and which method to use - Lax-Wendroff or RK4
    Limiter_1d(const double n_points, const double cfl,
               string method, const double running_time,
               string initial_data_indicator)
               :Finite_Volume_Solver_1d(n_points,cfl,
               method, running_time,initial_data_indicator){}
    void run(); //true when output is to be given, and false when it doesn't;
protected:
    void rhs_function();

    void foup(); //First order upwind scheme  
    void ssp_rk2_solver();
    void ssp_rk3_solver();

    void rhs_soup();
    void rhs_soup_minmod(); //This gives RHS of the system of ODEs and stores it to where rhs points.

    void soup_rk2(); //Second order upwind Scheme
    void soup_rk3();
    void soup_rk2_minmod(); //Second order upwind Scheme with min-mod limiter
    void soup_rk3_minmod(); //Second order upwind Scheme with min-mod limiter

    double limiter(int i); //psi(r_i)
};


void Limiter_1d::foup() //first order upwind scheme
{
  solution[0] = (1.0 - sigma) * solution_old[0] + cfl * solution_old[n_points - 1];
  for (int i = 1; i < n_points; i++)
  {
    solution[i] = (1.0 - sigma) * solution_old[i] + cfl * solution_old[i - 1];
  }
}

double Limiter_1d::limiter(int i) //Psi(r_i)
{

  double previous_value; //solution_old[i-1] isolated for handling i=0 case
  //while using the minimum lines
  if (i == 0)
    previous_value = solution_old[n_points-1];
  else
    previous_value = solution_old[i-1];
  if ((solution_old[i+1] - solution_old[i]) * (solution_old[i] -previous_value) <= 0.0) 
    return 0.0; //we are defining phi(a/b). If a or b is zero, phi is zero.
  else 
    {
    double r = (solution_old[i+1] - solution_old[i]) / (solution_old[i] -previous_value);
    if (method=="soup_rk3_minmod" || method=="soup_rk2_minmod")
      {
        return min(r,1.0);
      }
    else if (method =="soup_rk3_superbee" || method == "soup_rk2_superbee")
      {
        double value = max(0.0,min(2.0*r,1.0));
        value = max(value,min(r,2.0));
        return value;
      }
    else
      {
        cout << "Incorrect method specified" <<endl;
        assert(false);
      }
    }
}

void Limiter_1d::rhs_function()
{
  //rhs_function changes with the scheme.

  //first scheme is soup
  //rhs = -a/h(3/2 v_i - 2 v_{i-1} + 1/2 v_{i-2})
  if (method == "soup_rk2" || method == "soup_rk3")
  {
    /*
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
    ul = soup_reconstructor(solution[n_points-1], solution[n_points-2]);//u_{-1/2}^L
    flux = coefficient*ul; //f_{-1/2}
    (*rhs)[0]          +=  flux / h; //rhs[0] = -(f_{1/2} - f_{-1/2})/h,
    (*rhs)[n_points-1] += -flux / h;
    //rhs[n-1]=-(f_{n-1/2}-f_{n-3/2})=-(f_{-1/2}-f_{n-3/2})
    

    //Next compute f_{1/2}
    ul = soup_reconstructor(solution[0],solution[n_points-1]); //u_{1/2}^L
    flux = coefficient*ul; //f_{1/2}
    (*rhs)[1] +=  flux/h; //rhs[1] =-(f_{3/2}-f_{1/2})/h
    (*rhs)[0] += -flux/h; //rhs[0] =-(f_{1/2}-f_{-1/2})/h
    for (int i = 2; i<n_points;i++) //Only n_points-1 fluxes to be computed
    {
      ul = soup_reconstructor(solution[i-1], solution[i-2]); //u_{i-1/2}^L
      flux = coefficient * ul;//f_{i-1/2}
      (*rhs)[i]  +=  flux/h; //rhs[i]  =-(f_{i+1/2}-f_{i-1/2})/h
      (*rhs)[i-1]+= -flux/h; //rhs[i-1]=-(f_{i-1/2}-f_{i-3/2})/h
    }*/
    
    (*rhs)[0] = - (coefficient/h) * (1.5 * solution[0] 
                                 - 2.0 * solution[n_points-1] 
                                 + 0.5 * solution[n_points-2]);
    (*rhs)[1] = - (coefficient/h) * (1.5 * solution[1] 
                                 - 2.0 * solution[0] 
                                 + 0.5 * solution[n_points-1]);
    for (int i = 2; i < n_points; i++)
    {
        (*rhs)[i] = - (coefficient/h) * (1.5 * solution[i] 
                                    - 2.0 * solution[i-1] 
                                    + 0.5 * solution[i-2]);
    }
  }
  //second scheme is soup with minmod limiter
  //rhs = -a/h[v_i - v_{i-1} + 1/2( phi(r_i)(v_i-v_{i-1}) - phi(r_{i-1})(v_{i-1} - v_{i-2}) )]
  else if (method == "soup_rk3_minmod" || method == "soup_rk2_minmod" ||
           method == "soup_rk3_superbee" || method == "soup_rk2_superbee") //limiter will change according to method
  {
    (*rhs)[0] = -(coefficient/h) * ( (solution[0] - solution[n_points-1])
                                  +0.5  * limiter(0)  //Phi(r_0)
                                        * (solution[0] - solution[n_points-1])
                                  -0.5  * limiter(n_points-1)
                                        * (solution[n_points-1] - solution[n_points-2]));
    (*rhs)[1] = -(coefficient/h) * ( (solution[1] - solution[0])
                                  +0.5  * limiter(1)  //Phi(r_1)
                                        * (solution[1] - solution[0])
                                  -0.5  * limiter(0)
                                        * (solution[0] - solution[n_points-1]));
    for (int i = 2; i<n_points; i++)
    {
      (*rhs)[i] = -(coefficient/h) * ( (solution[i] - solution[i-1])
                                  +0.5 * limiter(i) //Phi(r_i)
                                        * (solution[i] - solution[i-1])
                                  -0.5 * limiter(i-1)
                                        * (solution[i-1] - solution[i-2]));
    }
  }
  else 
  {
    cout << "Incorrect method in rhs_function(method)" << endl;
    assert(false);
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
    /*
    Can't find the bug in this
    rhs = &temp;
    rhs_function(); //Put temp = rhs(solution)=rhs(solution_old)
    add(solution_old,dt,temp,k1); //Put k1 = solution_old + dt*rhs(solution_old)
    
    //We want to put k2 = 3/4 * U^n + 1/4*(k1 + dt*rhs(k1))
    solution = k1; //Since this is the thing that rhs() picks up
    rhs = &temp;
    rhs_function(); //Now, temp = rhs(solution) = rhs(k1)
    add(k1, dt, temp,solution); //Put solution = k1 + dt*rhs(k1)
    for (int i = 0; i<n_points; i++) 
    {
      solution[i] = 0.25 * solution[i];
    } //Now solution = 1/4(k1 + dt * rhs(k1))
    add(solution,0.75*dt,solution_old,k2);
    rhs = &temp;
    solution = k2;
    rhs_function(); //puts temp = rhs(solution)=rhs(k2)
    for (int i = 0; i<n_points; i++) 
    {
      solution[i] = (1.0/3.0) * solution_old[i] + (2.0/3.0)*(k2[i] + dt * temp[i]);
    }
    */
    /*
    rhs = &temp;
    //k1 = rhs_function(solution) = rhs_function(solution_old)
    rhs_function();//depends on method
    //Temporarily putting u^{n+1} = solution_new = solution_old + dt/2 * k1
    add(solution_old, dt / 3.0, temp, solution);
    //So, computing k2 = rhs_function(u^n + dt/2 * k1)
    rhs_function();
    add(solution_old, 0.5 * dt, temp, solution);
    rhs_function();
    add(solution_old, dt, temp, solution);
    */
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
