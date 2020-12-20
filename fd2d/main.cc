#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <string>
#include <stdio.h>
#include <sys/time.h>

#include "array2d.h"
#include "vtk_anim.h"
using namespace std;

//Returns true if real number is integer, false otherwise.
bool int_tester(double a)
{
    double dummy;
    double c = modf(a,&a);
    //cout << "Fractional part is "<<c<<endl;
    return (c<1e-4);
}

class Linear_Convection_2d
{
public:
    Linear_Convection_2d(const double n_points,
                         double lam_x, /*const double lam_y, */
                         string method, const double running_time,
                         int initial_data_indicator); 

    void run();
    void get_error(vector<double> &l1_vector, vector<double> &l2_vector, 
                   vector<double> &linfty_vector, vector<double> &snapshot_error);
private:
    void make_grid();
    void set_initial_solution();

    void upwind();
    void lw();
    void ct_upwind();

    double hat_function(double grid_point);
    double step_function(double grid_point); 
    double exp_func_25(double grid_point);

    void evaluate_error_and_output_solution(const int time_step_number);
    vector<double> grid_x,grid_y;
    double theta, coefficient_x, coefficient_y, x_min, x_max, y_min, y_max;

    double interval_part(double);

    Array2D solution_old; //Solution at previous step
    Array2D solution; //Solution at present step
    Array2D initial_solution;

    Array2D solution_exact; //Exact solution at present time step

    Array2D error;
    double snapshot_error;
    //After a certain time, by periodicity, the solution equals the initial soln
    //For the PDE qt + uqx + vqy = 0, the exact solution is q(x,y,t)=f(x-ut,y-vt)
    //Since we are assuming periodicity in both $x$ and $y$ directions with
    //period 2, solution agrees with initial solution when u*t, v*t are even integers.
    //So, we return to original solution for t*u = 2n, t*v= 2m for some n,m. 
    //We cannot always ensure this. If u = sqrt(2), v = sqrt(3), for u*t
    //to be an integer, we need t = m/sqrt(2), but then v*t cannot 
    //be an integer. Thus, the condition is, that u/v is a rational number

    double n_points, dx, dy, dt, t, running_time;
    double lam_x,lam_y, sigma_x,sigma_y;

    string method;
    int initial_data_indicator;
};

Linear_Convection_2d::Linear_Convection_2d(double n_points, 
                                           double lam_x,/* const double lam_y,*/
                                           string method,
                                           double running_time, int initial_data_indicator):
                                           n_points(n_points), 
                                           lam_x(lam_x),
                                           running_time(running_time), method(method),
                                           initial_data_indicator(initial_data_indicator)
{
    theta = M_PI/4.0;
    coefficient_x = 1.0, coefficient_y = 1.0;
    //coefficient_x = 1.0, coefficient_y = 1.0;
    x_min = -1.0, x_max = 1.0, y_min = -1.0, y_max = 1.0;
    dx = (x_max - x_min) / (n_points), dy = (y_max-y_min)/(n_points);
    //In interval [0,1], if we run the loop for i = 0,1,...,n-1
    //and take the grid spacing to be 1/n, we won't reach the end of interval.
    t = 0.0;
    dt = lam_x*dx/(abs(coefficient_x));
    /*//We want dt/coefficient_x to be integer, this ensures it.
    dt = dt/coefficient_x; Or maybe this step is not needed since dt already
    divided by coefficient_x
    */
    lam_x = abs(coefficient_x)*dt/dx;//lam_x updates
    sigma_x = coefficient_x*dt/(dx), sigma_y = coefficient_y*dt/(dy);
    //Since we have lam_x = |u|dt/dx,lam_y = |v|dt/dy, we can actually write
    //Thus, dt = lam_x*dx/|u|. Thus, lam_y = |v/u|*(dx/dy) * lam_x = 
    lam_y = (coefficient_y/coefficient_x)*dx/dy * lam_x;
    if (method == "lw"    &&
        (dx - dy < 1e-12) && 
        (dt/dx > 1.0/sqrt(coefficient_x*coefficient_x + coefficient_y*coefficient_y)))
    {
      cout << "WARNING - You are using lw with unstable sigma_x = "<< sigma_x <<endl;
    }
    snapshot_error = -1.0; //Put bad value for testing success.
    cout << "dx = " << dx << endl;
    cout << "dy = " << dy << endl;
    cout << "dt = " << dt << endl;
    cout << "lam_x = " <<lam_x << endl;
    cout << "lam_y = " <<lam_y << endl;
    //grid.resize(n_points);
    error.resize(n_points,n_points);
    grid_x.resize(n_points),grid_y.resize(n_points);
    initial_solution.resize(n_points,n_points);
    solution_old.resize(n_points,n_points);
    solution.resize(n_points,n_points);
    solution_exact.resize(n_points,n_points);
}

void Linear_Convection_2d::make_grid()
{
  //Note that you must run two for loops for a rectangular grid.
  for (unsigned int i = 0; i < n_points; i++)
  {
    grid_x[i] = x_min + i * dx;
    grid_y[i] = y_min + i * dy;
  }
}

void Linear_Convection_2d::set_initial_solution()
{
    for (unsigned int i = 0; i < n_points; i++)
      for (unsigned int j = 0; j < n_points; j++)
      {
          double x = x_min + i*dx;
          double y = y_min + j*dy;
          switch (initial_data_indicator)
          {
          case 0:
              solution(i,j) = sin(2.0 * M_PI * x / (x_max - x_min))
                              * sin(2.0 * M_PI * y / (y_max - y_min));
              break;
          case 1:
              solution(i,j) = hat_function(x)*hat_function(y);
              break;
          case 2:
              solution(i,j) = step_function(x)*step_function(y);
              break;
          case 3:
              solution(i,j) = exp_func_25(x)*exp_func_25(y);
              break;
          default:
              cout << "You entered the wrong initial_data_indicator ";
              assert(false);
          }
          initial_solution(i,j) = solution(i,j);
      }
}

void Linear_Convection_2d::upwind()
{
  //(i,j)=(0,0)
  solution(0,0) = (1.0-lam_x-lam_y)*solution_old(0,0)
                    +max(coefficient_x,0.)*(dt/dx)*solution_old(n_points-1,0)
                    +max(coefficient_y,0.)*(dt/dy)*solution_old(0,n_points-1)
                    -min(coefficient_x,0.)*(dt/dx)*solution_old(1,0)
                    -min(coefficient_y,0.)*(dt/dy)*solution_old(0,1);
  //i=0
  for (unsigned int j = 1; j<n_points-1;j++) 
    {
      solution(0,j) = (1.0-lam_x-lam_y)*solution_old(0,j)
                      +max(coefficient_x,0.)*(dt/dx)*solution_old(n_points-1,j)
                      +max(coefficient_y,0.)*(dt/dy)*solution_old(0,j-1)
                      -min(coefficient_x,0.)*(dt/dx)*solution_old(1,j)
                      -min(coefficient_y,0.)*(dt/dy)*solution_old(0,j+1);
    }
    //i = 0; j =n_points-1
  solution(0,n_points-1) = (1.0-lam_x-lam_y)*solution_old(0,n_points-1)
                          +max(coefficient_x,0.)*(dt/dx)*solution_old(n_points-1,n_points-1)
                          +max(coefficient_y,0.)*(dt/dy)*solution_old(0,n_points-2)
                          -min(coefficient_x,0.)*(dt/dx)*solution_old(1,n_points-1)
                          -min(coefficient_y,0.)*(dt/dy)*solution_old(0,0);
  //j=0
  for (unsigned int i = 1; i<n_points-1 ;i++) 
  {
    solution(i,0) = (1.0-lam_x-lam_y)*solution_old(i,0)
                        +max(coefficient_x,0.)*(dt/dx)*solution_old(i-1,0)
                        +max(coefficient_y,0.)*(dt/dy)*solution_old(i,n_points-1)
                        -min(coefficient_x,0.)*(dt/dx)*solution_old(i+1,0)
                        -min(coefficient_y,0.)*(dt/dy)*solution_old(i,1);;
  }
  //j=0, i = n_points-1
  solution(n_points-1,0) = (1.0-lam_x-lam_y)*solution_old(n_points-1,0)
                  +max(coefficient_x,0.)*(dt/dx)*solution_old(n_points-2,0)
                  +max(coefficient_y,0.)*(dt/dy)*solution_old(n_points-1,n_points-1)
                  -min(coefficient_x,0.)*(dt/dx)*solution_old(0,0)
                  -min(coefficient_y,0.)*(dt/dy)*solution_old(n_points-1,1);
  //i = n_points-1
  for (unsigned int j = 1;j<n_points-1;j++)
  {
    solution(n_points-1,j) = (1.0-lam_x-lam_y)*solution_old(n_points-1,j)
                          +max(coefficient_x,0.)*(dt/dx)*solution_old(n_points-2,j)
                          +max(coefficient_y,0.)*(dt/dy)*solution_old(n_points-1,j-1)
                          -min(coefficient_x,0.)*(dt/dx)*solution_old(0,j)
                          -min(coefficient_y,0.)*(dt/dy)*solution_old(n_points-1,j+1);
  }
  //j = n_points - 1
  for (unsigned int i = 1; i<n_points-1; i++)
  {
        solution(i,n_points-1) = (1.0-lam_x-lam_y)*solution_old(i,n_points-1)
                        +max(coefficient_x,0.)*(dt/dx)*solution_old(i-1,n_points-1)
                        +max(coefficient_y,0.)*(dt/dy)*solution_old(i,n_points-2)
                        -min(coefficient_x,0.)*(dt/dx)*solution_old(i+1,n_points-1)
                        -min(coefficient_y,0.)*(dt/dy)*solution_old(i,0);
  }
  //i = n_points - 1, j = n_points-1
  solution(n_points-1,n_points-1) = (1.0-lam_x-lam_y)*solution_old(n_points-1,n_points-1)
                +max(coefficient_x,0.)*(dt/dx)*solution_old(n_points-2,n_points-1)
                +max(coefficient_y,0.)*(dt/dy)*solution_old(n_points-1,n_points-2)
                -min(coefficient_x,0.)*(dt/dx)*solution_old(0,n_points-1)
                -min(coefficient_y,0.)*(dt/dy)*solution_old(n_points-1,0);
  //i != 0,n_points-1, j != 0,n_points-1
  for (unsigned int i = 1; i < n_points-1; i++)
      for (unsigned int j = 1; j < n_points-1; j++)
      {
        solution(i,j) = (1.0-lam_x-lam_y)*solution_old(i,j)
                        +max(coefficient_x,0.)*(dt/dx)*solution_old(i-1,j)
                        +max(coefficient_y,0.)*(dt/dy)*solution_old(i,j-1)
                        -min(coefficient_x,0.)*(dt/dx)*solution_old(i+1,j)
                        -min(coefficient_y,0.)*(dt/dy)*solution_old(i,j+1);
      }
}

void Linear_Convection_2d::ct_upwind()
{
  //(i,j)=(0,0)
  solution(0,0) = (1.0-sigma_x)*(1.0-sigma_y)*solution_old(0,0)
                +sigma_x*(1-sigma_y)*solution_old(n_points-1,0)
                +(1-sigma_x)*sigma_y*solution_old(0,n_points-1)
                +sigma_x*sigma_y*solution_old(n_points-1,n_points-1);
  //i=0
  for (unsigned int j = 1; j<n_points;j++) 
    {
    solution(0,j) = (1.0-sigma_x)*(1.0-sigma_y)*solution_old(0,j)
                    +sigma_x*(1-sigma_y)*solution_old(n_points-1,j)
                    +(1-sigma_x)*sigma_y*solution_old(0,j-1)
                    +sigma_x*sigma_y*solution_old(n_points-1,j-1);
    }
  //Extremely stupid mistake - I thought that the above cover all
  //the extra cases, when it clearly doesn't.
  //j=0
  for (unsigned int i = 1; i<n_points ;i++) 
  {
  solution(i,0) = (1.0-sigma_x)*(1.0-sigma_y)*solution_old(i,0)
                  +sigma_x*(1-sigma_y)*solution_old(i-1,0)
                  +(1-sigma_x)*sigma_y*solution_old(i,n_points-1)
                  +sigma_x*sigma_y*solution_old(i-1,n_points-1);
  }
  for (unsigned int i = 1; i < n_points; i++)
      for (unsigned int j = 1; j < n_points; j++)
      {
        solution(i,j) = (1.0-sigma_x)*(1.0-sigma_y)*solution_old(i,j)
                        +sigma_x*(1-sigma_y)*solution_old(i-1,j)
                        +(1-sigma_x)*sigma_y*solution_old(i,j-1)
                        +sigma_x*sigma_y*solution_old(i-1,j-1);
      }
}

void Linear_Convection_2d::lw()
{
  //(i,j)=(0,0)
  solution(0,0) = solution_old(0,0)
                      -0.5*sigma_x*(solution_old(1,0)-solution_old(n_points-1,0))
                      -0.5*sigma_y*(solution_old(0,1)-solution_old(0,n_points-1))
                      +0.5*sigma_x*sigma_x*(solution_old(1,0)
                                            -2.0*solution_old(0,0)
                                            +solution_old(n_points-1,0))
                      +0.25*sigma_x*sigma_y*(solution_old(1,1)
                                            -solution_old(1,n_points-1)
                                            -solution_old(n_points-1,1)
                                            +solution_old(n_points-1,n_points-1))
                      +0.5*sigma_y*sigma_y*(solution_old(0,1)
                                            -2.0*solution_old(0,0)
                                            +solution_old(0,n_points-1));
  //i=0
  for (unsigned int j = 1; j<n_points-1;j++) 
    {
      solution(0,j) = solution_old(0,j)
                      -0.5*sigma_x*(solution_old(1,j)-solution_old(n_points-1,j))
                      -0.5*sigma_y*(solution_old(0,j+1)-solution_old(0,j-1))
                      +0.5*sigma_x*sigma_x*(solution_old(1,j)
                                            -2.0*solution_old(0,j)
                                            +solution_old(n_points-1,j))
                      +0.25*sigma_x*sigma_y*(solution_old(1,j+1)
                                            -solution_old(1,j-1)
                                            -solution_old(n_points-1,j+1)
                                            +solution_old(n_points-1,j-1))
                      +0.5*sigma_y*sigma_y*(solution_old(0,j+1)
                                            -2.0*solution_old(0,j)
                                            +solution_old(0,j-1));
    }
  //Extremely stupid mistake - I thought that the above cover all
  //the extra cases, when it clearly doesn't.
  //j=0
  for (unsigned int i = 1; i<n_points-1;i++) 
  {
    solution(i,0) = solution_old(i,0)
                -0.5*sigma_x*(solution_old(i+1,0)-solution_old(i-1,0))
                -0.5*sigma_y*(solution_old(i,1)-solution_old(i,n_points-1))
                +0.5*sigma_x*sigma_x*(solution_old(i+1,0)
                                      -2.0*solution_old(i,0)
                                      +solution_old(i-1,0))
                +0.25*sigma_x*sigma_y*(solution_old(i+1,1)
                                      -solution_old(i+1,n_points-1)
                                      -solution_old(i-1,1)
                                      +solution_old(i-1,n_points-1))
                +0.5*sigma_y*sigma_y*(solution_old(i,1)
                                      -2.0*solution_old(i,0)
                                      +solution_old(i,n_points-1));
  }
  //i=n_points-1,j=0
  solution(n_points-1,0) = solution_old(n_points-1,0)
                      -0.5*sigma_x*(solution_old(0,0)
                                    -solution_old(n_points-2,0))
                      -0.5*sigma_y*(solution_old(n_points-1,1)
                                    - solution_old(n_points-1,n_points-1))
                      +0.5*sigma_x*sigma_x*(solution_old(0,0)
                                            -2.0*solution_old(n_points-1,0)
                                            +solution_old(n_points-2,0))
                      +0.25*sigma_x*sigma_y*(solution_old(0,1)
                                            -solution_old(0,n_points-1)
                                            -solution_old(n_points-2,1)
                                            +solution_old(n_points-2,n_points-1))
                      +0.5*sigma_y*sigma_y*(solution_old(n_points-1,1)
                                            -2.0*solution_old(n_points-1,0)
                                            +solution_old(n_points-1,n_points-1));
  //i=0,j=n_points-1
  solution(0,n_points -1) = solution_old(0,n_points -1 )
                    -0.5*sigma_x*(solution_old(1,n_points -1 )
                                  -solution_old(n_points-1,n_points -1 ))
                    -0.5*sigma_y*(solution_old(0,0)-solution_old(0,n_points -2))
                    +0.5*sigma_x*sigma_x*(solution_old(1,n_points -1 )
                                          -2.0*solution_old(0,n_points -1 )
                                          +solution_old(n_points-1,n_points -1 ))
                    +0.25*sigma_x*sigma_y*(solution_old(1,0)
                                          -solution_old(1,n_points -2)
                                          -solution_old(n_points-1,0)
                                          +solution_old(n_points-1,n_points -2))
                    +0.5*sigma_y*sigma_y*(solution_old(0,0)
                                          -2.0*solution_old(0,n_points -1)
                                          +solution_old(0,n_points -2)); 
  //i=n_points-1
  for (unsigned int j = 1; j<n_points-1; j++)
  {
    solution(n_points-1,j) = solution_old(n_points-1,j)
                        -0.5*sigma_x*(solution_old(0,j)
                                      -solution_old(n_points-2,j))
                        -0.5*sigma_y*(solution_old(n_points-1,j+1)
                                      - solution_old(n_points-1,j-1))
                        +0.5*sigma_x*sigma_x*(solution_old(0,j)
                                              -2.0*solution_old(n_points-1,j)
                                              +solution_old(n_points-2,j))
                        +0.25*sigma_x*sigma_y*(solution_old(0,j+1)
                                              -solution_old(0,j-1)
                                              -solution_old(n_points-2,j+1)
                                              +solution_old(n_points-2,j-1))
                        +0.5*sigma_y*sigma_y*(solution_old(n_points-1,j+1)
                                              -2.0*solution_old(n_points-1,j)
                                              +solution_old(n_points-1,j-1));
  }
  // j = n_points -1 
  for (unsigned int i = 1; i<n_points -1 ; i++)
  {
    solution(i,n_points -1) = solution_old(i,n_points -1 )
                    -0.5*sigma_x*(solution_old(i+1,n_points -1 )
                                  -solution_old(i-1,n_points -1 ))
                    -0.5*sigma_y*(solution_old(i,0)-solution_old(i,n_points -2))
                    +0.5*sigma_x*sigma_x*(solution_old(i+1,n_points -1 )
                                          -2.0*solution_old(i,n_points -1 )
                                          +solution_old(i-1,n_points -1 ))
                    +0.25*sigma_x*sigma_y*(solution_old(i+1,0)
                                          -solution_old(i+1,n_points -2)
                                          -solution_old(i-1,0)
                                          +solution_old(i-1,n_points -2))
                    +0.5*sigma_y*sigma_y*(solution_old(i,0)
                                          -2.0*solution_old(i,n_points -1)
                                          +solution_old(i,n_points -2)); 
  }
  //(i,j) = (n_points-1,n_points-1)
  solution(n_points -1,n_points -1) = solution_old(n_points -1,n_points -1 )
                  -0.5*sigma_x*(solution_old(0,n_points -1 )
                                -solution_old(n_points -2,n_points -1 ))
                  -0.5*sigma_y*(solution_old(n_points -1,0)
                                -solution_old(n_points -1,n_points -2))
                  +0.5*sigma_x*sigma_x*(solution_old(0,n_points -1 )
                                        -2.0*solution_old(n_points -1,n_points -1 )
                                        +solution_old(n_points -2,n_points -1 ))
                  +0.25*sigma_x*sigma_y*(solution_old(0,0)
                                        -solution_old(0,n_points -2)
                                        -solution_old(n_points -2,0)
                                        +solution_old(n_points -2,n_points -2))
                  +0.5*sigma_y*sigma_y*(solution_old(n_points -1,0)
                                        -2.0*solution_old(n_points -1,n_points -1)
                                        +solution_old(n_points -1,n_points -2)); 
  for (unsigned int i = 1; i < n_points-1; i++)
      for (unsigned int j = 1; j < n_points-1; j++)
      {
        solution(i,j) = solution_old(i,j)
                        -0.5*sigma_x*(solution_old(i+1,j)-solution_old(i-1,j))
                        -0.5*sigma_y*(solution_old(i,j+1)-solution_old(i,j-1))
                        +0.5*sigma_x*sigma_x*(solution_old(i+1,j)
                                              -2.0*solution_old(i,j)
                                              +solution_old(i-1,j))
                        +0.25*sigma_x*sigma_y*(solution_old(i+1,j+1)
                                              -solution_old(i+1,j-1)
                                              -solution_old(i-1,j+1)
                                              +solution_old(i-1,j-1))
                        +0.5*sigma_y*sigma_y*(solution_old(i,j+1)
                                              -2.0*solution_old(i,j)
                                              +solution_old(i,j-1));
      }
}

void Linear_Convection_2d::evaluate_error_and_output_solution(int time_step_number)
{
    for (unsigned int i = 0; i < n_points; i++)
      for (unsigned int j = 0; j < n_points; j++)
      {
          double x = x_min + i*dx;
          double y = y_min + j*dy;
          switch (initial_data_indicator)
          {
          case 0:
              solution_exact(i,j) = sin(2.0 * M_PI * (x - coefficient_x*t) / (x_max - x_min))
                                    * sin(2.0 * M_PI * (y - coefficient_y*t) / (y_max - y_min));
              break;
          case 1:
              solution_exact(i,j) = hat_function(x - coefficient_x*t)
                                    * hat_function(y - coefficient_y*t);
              break;
          case 2:
              solution_exact(i,j) = step_function(x - coefficient_x*t)
                                    * step_function(y - coefficient_y*t);
              break;
          case 3:
              solution_exact(i,j) = exp_func_25(x-coefficient_x*t)
                                    * exp_func_25(y-coefficient_y*t);
              break;
          default:
              cout << "You entered the wrong initial_data_indicator ";
              assert(false);
          }
      }
    vtk_anim_sol(grid_x,grid_y,
          solution,
          t, time_step_number,
          "approximate_solution");
    vtk_anim_sol(grid_x,grid_y,
          solution_exact,
          t, time_step_number,
          "exact_solution");
    for (unsigned int i = 0; i < n_points; i++)
      for (unsigned int j = 0; j < n_points; j++)
      {
          error(i,j) = abs(solution(i,j) - solution_exact(i,j));
          /*cout << "Error at x = "<< x_min + i*dx;
          cout << ", y = "<< y_min + j*dy<< " = "<< error(i,j) << endl;*/
      }
}

void Linear_Convection_2d::run()
{
    make_grid();
    set_initial_solution(); //sets solution to be the initial data
    int time_step_number = 0; 
    evaluate_error_and_output_solution(time_step_number);
    while (t < running_time) //compute solution at next time step using solution_old
    {
        
        solution_old = solution;//update solution_old to be used at next time step
        if (method == "upwind")
            upwind();
        else if (method == "lw")
            lw();
        else if (method == "ct_upwind")
            ct_upwind();
        else
            assert(false);
        //Compute snapshot error for integert/coefficient_x, t/coefficient_y
        //You must run the solver for longer time, or you'd get very less error.
        if ((t>0.0)&&(int_tester(t*coefficient_x)==true)&&(int_tester(t*coefficient_y)))
        {
          //cout << "We are evaluating snapshot error at t = "<< t<<endl;
          for (unsigned int i = 0; i < n_points; i++)
            for (unsigned int j = 0; j < n_points; j++)
            {
                snapshot_error = max(snapshot_error,
                                     abs(solution(i,j) - initial_solution(i,j)));
                /*cout << "Error at x = "<< x_min + i*dx;
                cout << ", y = "<< y_min + j*dy<< " = "<< error(i,j) << endl;*/
            }
        }
        time_step_number += 1;
        t = t + dt; 
        evaluate_error_and_output_solution(time_step_number);
    }
    cout << "For n_points = " << n_points<<", we took ";
    cout << time_step_number << " steps." << endl;
}

void run_and_get_output(double n_points, double cfl,
                        string method, double running_time,
                        int initial_data_indicator,
                        unsigned int max_refinements)
{
    ofstream error_vs_h;
    error_vs_h.open("error_vs_h.txt");
    vector<double> linfty_vector;
    vector<double> l2_vector;
    vector<double> l1_vector;
    vector<double> snapshot_vector;

    for (unsigned int refinement_level = 0; refinement_level <= max_refinements;
        refinement_level++)
    {
        Linear_Convection_2d solver(n_points, cfl, method, running_time,
                                    initial_data_indicator);
        struct timeval begin, end;
        gettimeofday(&begin, 0);
        solver.run();
        gettimeofday(&end, 0);
        long seconds = end.tv_sec - begin.tv_sec;
        long microseconds = end.tv_usec - begin.tv_usec;
        double elapsed = seconds + microseconds * 1e-6;
        cout << "Time taken by this iteration is " << elapsed << " seconds." << endl;
        solver.get_error(l1_vector,l2_vector,linfty_vector,snapshot_vector);//push_back resp. error.
        cout << "Snapshot error is "<< snapshot_vector[refinement_level] << endl;
        error_vs_h << 2.0/n_points << " " << linfty_vector[refinement_level] << "\n";
        n_points = 2.0 * n_points;
        if (refinement_level > 0) //Computing convergence rate.
        {
            cout << "Linfty convergence rate at refinement level ";
            cout << refinement_level << " is ";
            cout << abs(log(linfty_vector[refinement_level] 
                             / linfty_vector[refinement_level - 1])) / log(2.0);
            cout << endl;
            
            cout << "L2 convergence rate at refinement level ";
            cout << refinement_level<< " is " ;
            cout << abs(log(l2_vector[refinement_level] 
                             / l2_vector[refinement_level - 1])) / log(2.0);
            cout << endl;

            cout << "L1 convergence rate at refinement level ";
            cout << refinement_level << " is " ;
            cout << abs(log(l1_vector[refinement_level] 
                             / l1_vector[refinement_level - 1])) / log(2.0);
            cout << endl;

            cout << "Snapshot convergence rate at refinement level ";
            cout << refinement_level << " is " ;
            cout << abs(log(snapshot_vector[refinement_level] 
                             / snapshot_vector[refinement_level - 1])) / log(2.0);
            cout << endl;
        }
        error_vs_h.close();
    }
    if (snapshot_vector[linfty_vector.size()-1]-(-1.0) < 1e-12)
    {
      cout<<"WARNING -Initial state never reached, so snapshot error not calculated.";
      cout<<" Check if running_time is too less, or coefficient_x/coefficient_y";
      cout<<" is irrational " <<endl;
    }
    cout << "After " << max_refinements << " refinements, l_infty error = ";
    cout << linfty_vector[linfty_vector.size()-1] << endl;
    cout << "The L2 error is " << l2_vector[linfty_vector.size()-1] << endl;
    cout << "The L1 error is " << l1_vector[linfty_vector.size()-1] << endl;
    cout << "The L_infty snapshot error is " <<snapshot_vector[linfty_vector.size()-1]<<endl;
    cout << endl;         
}

double Linear_Convection_2d::interval_part(double x)
{
    if (x > x_max)
        return x - ceil((x - x_max) / (x_max - x_min)) * (x_max - x_min);
    else if (x < x_min)
        return x + ceil((x_min - x) / (x_max - x_min)) * (x_max - x_min);
    else
        return x;
}

double Linear_Convection_2d::hat_function(double grid_point)
{
    double value;
    grid_point = interval_part(grid_point);
    if (grid_point < -0.5 || grid_point > 0.5)
        value = 0.0;
    else if (grid_point <= 0.0)
        value = grid_point + 0.5;
    else if (grid_point >= 0.0)
        value = 0.5 - grid_point;
    return value;
}

double Linear_Convection_2d::step_function(double grid_point)
{
    double value;
    grid_point = interval_part(grid_point);
    if (grid_point < x_min + (x_max - x_min) / 4.0 || 
        grid_point > x_max - (x_max - x_min) / 4.0)
        value = 0.0;
    else
        value = 1.0;
    return value;
}

double Linear_Convection_2d::exp_func_25(double grid_point)
{
  grid_point = interval_part(grid_point);
  return exp(-25*grid_point*grid_point);
}

void Linear_Convection_2d::get_error(vector<double> &l1_vector,
                                     vector<double> &l2_vector,
                                     vector<double> &linfty_vector,
                                     vector<double> &snapshot_vector)
{
    double l1 = 0.0;
    double l2 = 0.0;
    double linfty = 0.0;
    for (int j = 0; j < n_points; j++) 
      for (int i = 0; i<n_points; i++)
      {
          l1 = l1 + error(i,j)  * dx * dy;           // L1 error
          l2 = l2 + error(i,j) * error(i,j)  * dx * dy;// L2 error
          linfty = max(linfty, error(i,j));            // L_infty error
      }
    l2 = sqrt(l2);
    l1_vector.push_back(l1);
    l2_vector.push_back(l2);
    linfty_vector.push_back(linfty);
    snapshot_vector.push_back(snapshot_error);
}

int main(int argc, char **argv)
{
    if (argc != 6)
    {
        cout << "Incorrect format, use" << endl;
        cout << "./main upwind/lw/ct_upwind(for the respective method)";
        cout << " sigma_x running_time initial_data_indicator " ;
        cout << "max_refinements" << endl;
        cout << "Choices for method"<<endl;
        cout << "upwind"<<endl<<"lw"<<endl<<"ct_upwind"<<endl;
        cout << "Choices for initial data "<<endl;
        cout << "0 - smooth_sine"<<endl<<"1 - hat"<<endl;
	cout << "2 - step"<<endl;
        cout << "3 - exp_func_25" <<endl;
        assert(false);
    }
    string method = argv[1];
    cout << "method = " << method << endl;
    double n_points = 100.0;
    double cfl = stod(argv[2]);
    cout << "cfl = " << cfl << endl;
    double running_time = stod(argv[3]);
    cout << "running_time = " << running_time << endl;
    int initial_data_indicator = stoi(argv[4]);
    cout << "initial_data_indicator = " << initial_data_indicator << endl;
    unsigned int max_refinements = stoi(argv[5]);
    cout << "max_refinements = " << max_refinements <<endl;
    run_and_get_output(n_points, cfl, method, running_time,
                       initial_data_indicator, max_refinements);
}
