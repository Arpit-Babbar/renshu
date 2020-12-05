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

//Recall that f_{j+1/2} = a*u_{j+1/2}^L
//We roughly have u_{j+1/2}^L = u_j + phi(u_{j-1}-u_j,u_{j+1}-u_j)
//Where phi is a limiter like minmod, superbee, minmod3

//First we define the limiters phi.

//Roughly, this limiter gives min(abs(back_diff),abs(fwd_diff))
double minmod(double fwd_diff,double back_diff)
{
  if (fwd_diff * back_diff <= 0.0) 
      return 0.0;
  else 
      return min(fwd_diff / back_diff,1.0) * back_diff;
}

double superbee(double fwd_diff,double back_diff)
{
  if (fwd_diff * back_diff <= 0.0) 
      return 0.0;
  else 
      {
      double r = fwd_diff/back_diff;
      double temp = max(0.0,min(2.0*r,1.0));
      return max(temp, min(r,2.0)) * back_diff;
      }
}

//back_diff will mean uj - ujm1 for some j. Similar for cent_diff,fwd_diff.
double minmod3(double back_diff,double cent_diff,double fwd_diff)
{
  if ( (back_diff*cent_diff <= 0.0) ||
        (cent_diff*fwd_diff <= 0.0) )
      return 0.0;
  else //s min(|a|,|b|,|c|) where s = sign(a)=sign(b)=sign(c)
      {
        double temp = min(abs(back_diff),abs(cent_diff));
        temp = min(temp,abs(fwd_diff));
        return (back_diff/abs(back_diff)) * temp;
      }
}

//Here, we define the function that uses limiters to compute u_{j-1/2}^L
//where u_{j-1/2}^L = u_j + phi(...)
double reconstructor(double ujm2, double ujm1, double uj, string limiter)
{    
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

class Solver : public Finite_Volume_Solver_1d
{
public:
    //and which scheme to use - Lax-Wendroff or RK4
    Solver(const double n_points, const double cfl,
               string scheme, const double running_time,
               string initial_data_indicator, string limiter);
    void run(); //true when output is to be given, and false when it doesn't;
protected:
    void make_grid();
    void rhs_function();

    void foup(); //First order upwind scheme  
    void ssp_rk2_solver();
    void ssp_rk3_solver();

    void rhs_soup();

    string limiter;

    void soup_rk2(); //Second order upwind Scheme
    void soup_rk3();
};

Solver::Solver(const double n_points, const double cfl,
               string scheme, const double running_time,
               string initial_data_indicator, string limiter)
               :Finite_Volume_Solver_1d(n_points,cfl,
               scheme, running_time,initial_data_indicator), limiter(limiter){}

void Solver::make_grid()
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

void Solver::foup() //first order upwind scheme
{
  solution[0] = (1.0 - sigma) * solution_old[0] + cfl * solution_old[n_points - 1];
  for (int i = 1; i < n_points; i++)
  {
    solution[i] = (1.0 - sigma) * solution_old[i] + cfl * solution_old[i - 1];
  }
}

void Solver::rhs_function()
{
  fill((*rhs).begin(),(*rhs).end(),0.0); //Sets (*rhs) vector to zero.
  //There are n_points+1 points on which the flux needs to be evaluated
  //But, the last and first flux are the same, so we'd not evaluate them

  //These values of solution would be useful.
  
  //Also, each flux finds its use in computing the RHS at two of the 
  //'FVM grid' points, so we'd put it in both the places in same iteration.
  double flux;
  //Recall that we had defined flux_{j+1/2} = coefficient*u_{j+1/2}^L.
  //We'd just compute that u_{j+1/2}^L, and store it as ul.
  double ul;
  //For usual reasons, we need to consider the first and last flux separately
  //First we compute f_{-1/2}
  ul = reconstructor(solution[n_points-2],
                     solution[n_points-1],
                     solution[0],
                     limiter);//u_{0-1/2}^L = u_{-1/2}^L
  flux = coefficient*ul; //f_{-1/2}
  (*rhs)[0]          +=  flux / h; //rhs[0] = -(f_{1/2} - f_{-1/2})/h,
  (*rhs)[n_points-1] += -flux / h;
  //rhs[n-1]=-(f_{n-1/2}-f_{n-3/2})=-(f_{-1/2}-f_{n-3/2})
  
  //Next compute f_{1/2}
  ul = reconstructor(solution[n_points-1], 
                     solution[0],
                     solution[1],
                     limiter); //u_{1-1/2}^L=u_{1/2}^L
  flux = coefficient*ul; //f_{1/2}
  (*rhs)[1] +=  flux/h; //rhs[1] =-(f_{3/2}-f_{1/2})/h
  (*rhs)[0] += -flux/h; //rhs[0] =-(f_{1/2}-f_{-1/2})/h*/
  for (int j = 2; j<n_points;j++) //Only n_points-1 fluxes to be computed
  {
    ul = reconstructor(solution[j-2],solution[j-1],solution[j],
                       limiter); //u_{i-1/2}^L
    flux = coefficient * ul;//f_{i-1/2}
    (*rhs)[j]  +=  flux/h; //rhs[i]  =-(f_{i+1/2}-f_{i-1/2})/h
    (*rhs)[j-1]+= -flux/h; //rhs[i-1]=-(f_{i-1/2}-f_{i-3/2})/h
  }
}

void Solver::ssp_rk3_solver()
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

void Solver::ssp_rk2_solver()
{
    rhs = &temp;
    rhs_function();
    add(solution_old,dt,temp,solution);
    rhs_function();
    add(solution,dt,temp,temp);
    add(0.5,solution_old,0.5,temp,solution);    
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
      if (scheme == "lw")
        lax_wendroff();
      else if (scheme == "foup")
        foup();
      else if (scheme == "soup_rk3")
        ssp_rk3_solver();
      else if (scheme == "soup_rk2")
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

template<typename Pde_Solver>
void run_and_get_output(double n_points, double max_refinements,
                        Pde_Solver &solver) //Never forget to put & here!!
{
    ofstream error_vs_h;
    error_vs_h.open("error_vs_h.txt");
    vector<double> linfty_vector;
    vector<double> l2_vector;
    vector<double> l1_vector;

    for (unsigned int refinement_level = 0; refinement_level <= max_refinements;
        refinement_level++)
    {
        if (refinement_level>0)
          solver.refine();
        struct timeval begin, end;
        gettimeofday(&begin, 0);
        solver.run();
        gettimeofday(&end, 0);
        long seconds = end.tv_sec - begin.tv_sec;
        long microseconds = end.tv_usec - begin.tv_usec;
        double elapsed = seconds + microseconds * 1e-6;
        cout << "Time taken by this iteration is " << elapsed << " seconds." << endl;
        solver.get_error(l1_vector,l2_vector,linfty_vector);//push_back resp. error.
        error_vs_h << 2.0/n_points << " " << linfty_vector[refinement_level] << "\n";
        n_points = 2.0 * n_points;
        if (refinement_level > 0)
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
        }
        if (refinement_level == max_refinements + 1) solver.output_final_error();
    }
    error_vs_h.close();
    cout << "After " << max_refinements << " refinements, l_infty error = ";
    cout << linfty_vector[max_refinements-1] << endl;
    cout << "The L2 error is " << l2_vector[max_refinements-1] << endl;
    cout << "The L1 error is " << l1_vector[max_refinements-1] << endl;
    cout << endl;         
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
                                  initial_data_indicator, limiter);
  run_and_get_output(n_points, max_refinements, solver);
}

