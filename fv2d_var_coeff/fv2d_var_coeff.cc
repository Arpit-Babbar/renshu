//Solving Q_t + uQ_x + vQ_y = 0

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <string>
#include <stdio.h>
#include <sys/time.h>

#include "/mnt/c/Users/arpit/Documents/GitHub/arpit_practise/include/array2d.h"
#include "/mnt/c/Users/arpit/Documents/GitHub/arpit_practise/include/vtk_anim.h"
#include "/mnt/c/Users/arpit/Documents/GitHub/arpit_practise/include/initial_conditions.h"
using namespace std;

//Returns true if real number is integer, false otherwise.
bool int_tester(double a)
{
  double dummy;
  double c = modf(a,&a);
  //cout << "Fractional part is "<<c<<endl;
  return (c<1e-4);
}

double reconstruct(double Qm1,double Q0)
{
  return Q0;
}


void (*update_flux)(int n_x, int n_y,double vel[2],double Q_l,double Q_r, double &flux);

void upwind(int n_x, int n_y, double vel[2], double Q_l, double Q_r, double& flux)
{
  const double v_n = vel[0]*n_x + vel[1]*n_y;//normal velocity
  flux = max(vel[n_x,n_y],0.0)*Q_l+ min(vel[n_x,n_y],0.)*Q_r;
}

//Computes advection_velocity in x direction at (x,y)
void rotational_velocity(double x, double y, double vel[2])
{
  vel[0] = -y, vel[1] = x;
}

void constant_velocity(double x, double y, double vel[2])
{
  vel[0] = 1.0, vel[1] = 1.0;
}

void (*advection_velocity)(double, double, double vel[2]) = &rotational_velocity;

class Linear_Convection_2d
{
public:
    Linear_Convection_2d(double N_x, double N_y,
                         double lam,
                         string method, const double running_time,
                         int initial_data_indicator); 

    void run(bool output_indicator);
    void get_error(vector<double> &l1_vector, vector<double> &l2_vector, 
                   vector<double> &linfty_vector, vector<double> &snapshot_error);
private:
    void make_grid();
    void set_initial_solution();

    void compute_time_step();//This computes the time step dt.

    //Computes flux_x(i+1/2,j), flux_y(i,j+1/2)
    void lw_normal_flux();
    void lw_flux(int i, int j,int n_x, int n_y, double& flux);
    void lw(int i, int j, int n_x, int n_y, double& flux);
    void upwind(int i, int j, int n_x, int n_y, double& flux);

    void apply_fvm();
    void apply_lw();

    void evaluate_error_and_output_solution(const int time_step_number,
                                            bool output_indicator);
    vector<double> grid_x,grid_y;
    
    double vel[2]; //advection velocity vector
    double theta, xmin, xmax, ymin, ymax;

    Array2D solution_old; //Solution at previous step
    Array2D solution; //Solution at present step
    Array2D initial_solution;

    Array2D residual;

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

    double N_x,N_y, dx, dy, dt, t, running_time;
    double lam, sigma_x, sigma_y;
    void update_ghost_values();
    void use_ghost_values();
    string method;
    int initial_data_indicator;
    I_Functions initial_function;
};

Linear_Convection_2d::Linear_Convection_2d(double N_x, double N_y, 
                                           double lam,
                                           string method,
                                           double running_time, 
                                           int initial_data_indicator):
                                           N_x(N_x), N_y(N_y), 
                                           lam(lam),
                                           running_time(running_time),
                                           method(method),
                                           initial_data_indicator(initial_data_indicator)
{
    theta = M_PI/4.0;
    xmin = -1.0, xmax = 1.0, ymin = -1.0, ymax = 1.0;
    dx = (xmax - xmin) / (N_x), dy = (ymax-ymin)/(N_y);
    //In interval [0,1], if we run the loop for i = 0,1,...,n-1
    //and take the grid spacing to be 1/n, we won't reach the end of interval.
    t = 0.0;
    
    initial_function.set(initial_data_indicator,xmin,xmax,ymin,ymax);
    snapshot_error = -1.0; //Put bad value for testing success.
    cout << "dx = " << dx << endl;
    cout << "dy = " << dy << endl;
    cout << "lam = " <<lam << endl;

    error.resize(N_x,N_y);
    grid_x.resize(N_x),grid_y.resize(N_y);
    initial_solution.resize(N_x,N_y);
    solution_old.resize(N_x,N_y,1);
    solution.resize(N_x,N_y,1);
    residual.resize(N_x,N_y,1);
    solution_exact.resize(N_x,N_y);

}

void Linear_Convection_2d::compute_time_step()
{
  double u0max=0.,u1max=0.;
  double vel[2];
  for (int i = 0; i<N_x; i++)
    for (int j = 0; j<N_y;j++) //Loop over all cell centers.
    {
      double x = xmin + 0.5*dx + i*dx, y = ymin + 0.5*dy + j*dy;
      (*advection_velocity)(x,y,vel);
      u0max = max(u0max,abs(vel[0])), u1max = max(u1max,abs(vel[1]));
    }
  double c;
  if (method == "upwind")
    c = 1.0;
  else if (method == "lw")
    c = 0.72;
  else 
  {
    cout << "Incorrect method"<<endl;
    assert(false);
  }
  dt = lam*c/(u0max/dx+u1max/dy);
  cout << "dt = "<<dt <<endl;
}

//Computes flux at at the face with centre
// (x_{i+0.5*n_x},y_{j+0.5*n_y}). In this code, we'd just have
//(n_x,n_y) = (1,0) or (0,1), so the centres will just be 
//(x_{i+1/2,j},y_j) or (x_i,y_{j+1/2})
void Linear_Convection_2d::upwind(int i, int j,
                                  int n_x,int n_y,
                                  double& flux)
{
  flux = max(vel[n_x,n_y],0.)*solution_old(i,j)
           + min(vel[n_x,n_y],0.)*solution_old(i+n_x,j+n_y);
}

//This function does the actual job of computing the flux.
void Linear_Convection_2d::lw(int i, int j, int n_x, int n_y,
                              double& flux)
{
  flux = 0.5*(vel[n_x,n_y])*(solution_old(i,j) + solution_old(i+n_x,j+n_y))
        -0.5*(vel[n_x,n_y])*(vel[n_x,n_y])*(n_x*dt/dx+n_y*dt/dy)
            *(solution_old(i+n_x,j+n_y)- solution_old(i,j))
        -0.125*vel[n_x,n_y]*vel[n_y,n_x]*(n_x*dt/dy + n_y*dt/dx)
              *(solution_old(i+n_y,j+n_x)-solution_old(i-n_y,j-n_x)
                +solution_old(i+1,j+1)-solution_old(i+n_x-n_y,j-n_x+n_y));
}



void Linear_Convection_2d::make_grid()
{
  //Note that you must run two for loops for a rectangular grid.
  for (unsigned int i = 0; i < N_x; i++)
    grid_x[i] = (xmin+0.5*dx) + i * dx;
  for (unsigned int j = 0; j<N_y; j++)
    grid_y[j] = (ymin+0.5*dy) + j * dy;
}

void Linear_Convection_2d::update_ghost_values()
{
  //For readability, we do corners separately.
  solution_old(-1,-1)= solution_old(N_x-1,N_y-1);
  solution_old(N_x,N_y)  = solution_old(0,0);
  solution_old(N_x,-1) = solution_old(0,N_y-1);
  solution_old(-1,N_y) = solution_old(N_x-1,0);

  solution(-1,-1)= solution(N_x-1,N_y-1);
  solution(N_x,N_y)  = solution(0,0);
  solution(N_x,-1) = solution(0,N_y-1);
  solution(-1,N_y) = solution(N_x-1,0);
  for (int j = 0; j<N_y;j++) //not doing corners
  {
    solution_old(-1,j) = solution_old(N_x-1,j);
    solution_old(N_x,j)  = solution_old(0,j);
    solution(-1,j)   = solution(N_x-1,j);
    solution(N_x,j)  = solution(0,j);
  }
  //j = -1, N
  for (int i = 0; i<N_x;i++)  //not doing corners
  {
    solution_old(i,-1) = solution_old(i,N_y-1);
    solution_old(i,N_y)  = solution_old(i,0);
    solution(i,-1) = solution(i,N_y-1);
    solution(i,N_y)  = solution(i,0);
  }
}

void Linear_Convection_2d::set_initial_solution()
{
  double x,y;
  for (unsigned int i = 0; i < N_x; i++)
    for (unsigned int j = 0; j < N_y; j++)
    {
      x = (xmin+0.5*dx) + i*dx, y = (ymin+0.5*dy) + j*dy;
      solution(i,j) = initial_function.value(x,y);
      initial_solution(i,j) = solution(i,j); //Stored only
      //for snapshot error
    }
}

//dy/dt = res(u) 
//where the residual coming from FVM is given by
//res(u)_{i,j}^{n}=-(flux_x(i+1/2,j)-flux_x(i-1/2,j))/dx
//                 -(flux_y(i,j+1/2)-flux_y(i,j-1/2))/dy
//res(u)_{i+1,j}^{n}=-(flux_x(i+3/2,j)-flux_x(i+1/2,j))/dx
//                   -(flux_y(i+1,j+1/2)-flux_y(i+1,j-1/2))/dy
//res(u)_{i,j+1}^{n}=-(flux_x(i+1/2,j+1)-flux_x(i-1/2,j+1))/dx
//                   -(flux_y(i,j+3/2)-flux_y(i,j+1/2))/dy
//res(u)_{i+1,j+1}^{n}=-(flux_x(i+3/2,j+1)-flux_x(i+1/2,j+1))/dx
//                     -(flux_y(i+1,j+3/2) - flux_y(i+1,j+1/2))/dx
//We will run loops over faces centred at (x_{i+1/2},y_j) and (x_i,y_{j+1/2})
//to compute flux_x(i+1/2,j), flux_y(i,j+1/2) and move it to the respective
//values of the residual.
//We avoid repetetion of faces in two ways - 1)one flux shows up in two places
//with opposite signs(look at values above), we only compute it once and move it to 
//the correct place 2) Last and 0th flux are the same, that flux shows up twice
//with opposite signs and we handle it accordingly.

void Linear_Convection_2d::apply_fvm()
{
  //double x,y;
  double flux; //flux_x(i+1/2,j), flux_y(i,j+1/2)
  //This loop computes the fluxes and adds them to where they are needed
  update_ghost_values();
  residual = 0.0;//For different time integration
  lam = dt/(dx*dy);
  //We'd do solution = solution_old - dt/dx * (f_x(i+1/2,j)-f_x(i-1/2,j))
  //                                - dt/dx * (f_y(i,j+1/2)-f_y(i,j-1/2))

  //dy/dt = res(u)

  //flux in x direction computed and used to update solution
  //Basically, flux_x(i+1/2,j)

  //
  for (int i = 0; i < N_x; i++)
    for (int j = 0; j< N_y; j++)
    {
      double x = (xmin+dx)+i*dx, y = ymin+0.5*dy+j*dy; //Values on face centre
      //(x_{i+1/2},y_j)
      (*advection_velocity)(x,y,vel);
      const double Q_l = reconstruct(solution(i-1,j),solution(i,j));
      const double Q_r = reconstruct(solution(i,j),solution(i+1,j));
      (*update_flux)(1,0,vel,Q_l,Q_r,flux);
      residual(i,j)     += -flux*dy;
      if (i==N_x-1)
        residual(0,j)   +=  flux*dy;
      else
        residual(i+1,j) +=  flux*dy;
    }

  //flux in y direction computed and used to update solution
  //Basically, flux_y(i,j+1/2)
  for (int i = 0; i < N_x; i++)
    for (int j = 0;j< N_y; j++)
    {
      double x = (xmin+0.5*dx)+i*dx, y = (ymin+dy)+j*dy; //Values on face centre.
      //(x_i,y_{j+1/2})
      (*advection_velocity)(x,y,vel);
      const double Q_l = reconstruct(solution(i-1,j),solution(i,j));
      const double Q_r = reconstruct(solution(i,j),solution(i,j+1));
      (*update_flux)(0,1,vel,Q_l,Q_r,flux);
      residual(i,j)     += -flux*dx;
      if (j==N_y-1)
        residual(i,0)   +=  flux*dx;
      else
        residual(i,j+1) +=  flux*dx;
    }

  for (int i = 0; i<N_x; i++)
    for (int j = 0; j<N_y; j++)
    {
      solution(i,j) = solution_old(i,j) + lam*residual(i,j);
    }
}

void Linear_Convection_2d::apply_lw()
{
  //double x,y;
  double flux; //flux_x(i+1/2,j), flux_y(i,j+1/2)
  //This loop computes the fluxes and adds them to where they are needed
  update_ghost_values();
  residual = 0.0;//For different time integration
  lam = dt/(dx*dy);
  //We'd do solution = solution_old - dt/dx * (f_x(i+1/2,j)-f_x(i-1/2,j))
  //                                - dt/dx * (f_y(i,j+1/2)-f_y(i,j-1/2))

  //dy/dt = res(u)

  //flux in x direction computed and used to update solution
  //Basically, flux_x(i+1/2,j)

  //
  for (int i = 0; i < N_x; i++)
    for (int j = 0; j< N_y; j++)
    {
      double x = (xmin+dx)+i*dx, y = ymin+0.5*dy+j*dy; //Values on face centre
      //(x_{i+1/2},y_j)
      (*advection_velocity)(x,y,vel);
      lw(i,j,1,0,flux); //flux_x(i+1/2,j)
      residual(i,j)     += -flux*dy;
      if (i==N_x-1)
        residual(0,j)   +=  flux*dy;
      else
        residual(i+1,j) +=  flux*dy;
    }

  //flux in y direction computed and used to update solution
  //Basically, flux_y(i,j+1/2)
  for (int i = 0; i < N_x; i++)
    for (int j = 0;j< N_y; j++)
    {
      double x = (xmin+0.5*dx)+i*dx, y = (ymin+dy)+j*dy; //Values on face centre.
      //(x_i,y_{j+1/2})
      (*advection_velocity)(x,y,vel);
      lw(i,j,0,1,flux);//flux_y(i,j+1/2)
      //Recall that, in FVM, solution updates as
      //Q_{i,j}^{n+1} = Q_{i,j}^n-(flux_x(i+1/2,j)-flux_x(i-1/2,j))*(dt/dx)
      //                         -(flux_y(i,j+1/2)-flux_y(i,j-1/2))*(dt/dy)
      //Q_{i+1,j}^{n+1} = Q_{i+1,j}^n-(flux_x(i+3/2,j)-flux_x(i+1/2,j))*(dt/dx)
      //                         -(flux_y(i+1,j+1/2)-flux_y(i+1,j-1/2))*(dt/dy)
      //Q_{i,j+1}^{n+1} = Q_{i,j+1}^n-(flux_x(i+1/2,j+1)-flux_x(i-1/2,j+1))*(dt/dx)
      //                         -(flux_y(i,j+3/2)-flux_y(i,j+1/2))*(dt/dy)
      //Q_{i+1,j+1}^{n+1} = Q_{i+1,j+1}^n-(flux_x(i+3/2,j+1)-flux_x(i+1/2,j+1))*(dt/dx)
      //                         -(flux_y(i+1,j+3/2) - flux_y(i+1,j+1/2))*(dt/dx)
      residual(i,j)     += -flux*dx;
      if (j==N_y-1)
        residual(i,0)   +=  flux*dx;
      else
        residual(i,j+1) +=  flux*dx;
    }

  for (int i = 0; i<N_x; i++)
    for (int j = 0; j<N_y; j++)
    {
      solution(i,j) = solution_old(i,j) + lam*residual(i,j);
    }
}


void Linear_Convection_2d::evaluate_error_and_output_solution(int time_step_number,bool output_indicator)
{
  double x=0.0,y=0.0;
  bool constant_indicator=false;//This indicates whether the coefficients are 
  //constant or not, to help compute exact solution. We set it to false by 
  //default, it'd be updated to true if next condition is true
  advection_velocity(x,y,vel);
  if (vel[0]==1.0&&vel[1]==1.0)
    constant_indicator = true;
  for (unsigned int i = 0; i < N_x; i++)
    for (unsigned int j = 0; j < N_y; j++)
    {
      x = (xmin+0.5*dx) + i*dx, y = (ymin+0.5*dy) + j*dy;
      advection_velocity(x,y,vel);
      solution_exact(i,j) = initial_function.exact_value(x,y,t,vel,constant_indicator);
    }
  if (output_indicator==true && time_step_number%5==0)
  {
  vtk_anim_sol(grid_x,grid_y,
        solution, solution_exact,
        t, time_step_number,
        "approximate_solution");
  }
  for (unsigned int i = 0; i < N_x; i++)
    for (unsigned int j = 0; j < N_y; j++)
    {
      error(i,j) = abs(solution(i,j) - solution_exact(i,j));
    }
}

void Linear_Convection_2d::run(bool output_indicator)
{
  make_grid();
  int time_step_number = 0;
  compute_time_step();//Computes dt
  set_initial_solution(); //sets solution to be the initial data
  evaluate_error_and_output_solution(time_step_number,output_indicator);
  while (t < running_time) //compute solution at next time step using solution_old
  {
    solution_old = solution;//update solution_old for next time_step
    if (method == "lw")
      apply_lw();
    else 
      apply_fvm();
    time_step_number += 1;
    t = t + dt; 
    evaluate_error_and_output_solution(time_step_number, output_indicator);
  }
  cout << "For N_x = " << N_x<<", N_y = "<<N_y<<", we took ";
  cout << time_step_number << " steps." << endl;
  if (output_indicator)
    cout <<"We produce output in this refinement level\n";
}

void run_and_output(double N_x, double N_y, double cfl,
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
    Linear_Convection_2d solver(N_x, N_y, cfl, method, running_time,
                                initial_data_indicator);
    //We calculate time takenṣ in our refinement.
    struct timeval begin, end; 
    gettimeofday(&begin, 0);
    solver.run(refinement_level==max_refinements);//Output only last soln
    gettimeofday(&end, 0);
    long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
    double elapsed = seconds + microseconds * 1e-6;
    cout << "Time taken by this refinement level is " << elapsed << " seconds." << endl;
    solver.get_error(l1_vector,l2_vector,linfty_vector,snapshot_vector);//push_back resp. error.
    double h = 2.*sqrt(1./(N_x*N_x) +1./(N_y*N_y));
    error_vs_h << h << " " << linfty_vector[refinement_level] << "\n";
    N_x = 2.0 * N_x,N_y = 2.0*N_y;
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
      //Temporary solution to snapshot rate, would be fine once we
      //move run_and_output to solver class.
      if (snapshot_vector.size()>0)
      {
        cout << "Snapshot convergence rate at refinement level ";
        cout << refinement_level << " is " ;
        cout << abs(log(snapshot_vector[refinement_level] 
                        / snapshot_vector[refinement_level - 1])) / log(2.0);
        cout << endl;
      }
    }
    error_vs_h.close();
  }
  cout << "After " << max_refinements << " refinements, l_infty error = ";
  cout << linfty_vector[linfty_vector.size()-1] << endl;
  cout << "The L2 error is " << l2_vector[linfty_vector.size()-1] << endl;
  cout << "The L1 error is " << l1_vector[linfty_vector.size()-1] << endl;
  if (snapshot_vector.size()>0)
    cout << "The L_infty snapshot error is " <<snapshot_vector[linfty_vector.size()-1]<<endl;
  cout << endl;         
}

void Linear_Convection_2d::get_error(vector<double> &l1_vector,
                                     vector<double> &l2_vector,
                                     vector<double> &linfty_vector,
                                     vector<double> &snapshot_vector)
{
    double l1 = 0.0;
    double l2 = 0.0;
    double linfty = 0.0;
    for (int j = 0; j < N_y; j++) 
      for (int i = 0; i < N_x; i++)
      {
        l1 = l1 + error(i,j)  * dx * dy;           // L1 error
        l2 = l2 + error(i,j) * error(i,j)  * dx * dy;// L2 error
        linfty = max(linfty, error(i,j));            // L_infty error
      }
    l2 = sqrt(l2);
    l1_vector.push_back(l1);
    l2_vector.push_back(l2);
    linfty_vector.push_back(linfty);
    //Compute snapshot error if last final time t has even integer 
    // t/u, t/v
    //You must run the solver for longer time, or you'd get very less error.
    //DOESN'T WORK FOR THIS, NEED 2pi period logic here
    if (int_tester(running_time/(2.0*M_PI)) == true &&
        int_tester(running_time/(2.0*M_PI)) == true)
    {
      cout << "Snapshot error successfuly computed.\n";
      for (unsigned int i = 0; i < N_x; i++)
        for (unsigned int j = 0; j < N_y; j++)
        {
          snapshot_error = abs(solution(i,j) - initial_solution(i,j));
        }
      snapshot_vector.push_back(snapshot_error);
    }
    else 
    {
      {
        cout<<"Initial state not reached on final time, so snapshot error";
        cout<<" not calculated. \n";
      }
    }
}

int main(int argc, char **argv)
{
    if (argc != 6 && argc != 7)
    {
      cout << "Incorrect format, use" << endl;
      cout << "./fd2d method";
      cout << " sigma_x running_time initial_data_indicator " ;
      cout << "max_refinements" << endl;
      cout << "Choices for method"<<endl;
      cout << "upwind, lw, ct_upwind, m_roe"<<endl;
      cout << "Choices for initial data "<<endl;
      cout << "0 - smooth_sine"<<endl<<"1 - hat"<<endl;
      cout << "2 - step"<<endl;
      cout << "3 - exp_func_25" <<endl;
      cout << "4 - exp_func_100\n";
      cout << "5 - cts_sine\n";
      cout << "You can add a 'constant' at the end of above to test";
      cout << "constant coefficients case. \n";
      assert(false);
    }
    if (argc == 7)
    {
      if (string(argv[6]) != "constant")
      {
        cout <<"Last argument must be constant.\n";
        cout << "You put "<< argv[6] <<endl;
        assert(false);
      }
      else
      {
        advection_velocity = &constant_velocity;
        cout <<"Scheme will be run with constant (u,v)=(1,1)"<<endl;
      }
    }
    string method = argv[1];
    cout << "method = " << method << endl;
    int N_x = 75, N_y = 75;
    double sigma_x = stod(argv[2]);
    cout << "sigma_x = " << sigma_x << endl;
    double running_time = stod(argv[3]);
    cout << "running_time = " << running_time << endl;
    int initial_data_indicator = stoi(argv[4]);
    cout << "initial_data_indicator = " << initial_data_indicator << endl;
    signed int max_refinements = stoi(argv[5]);
    cout << "max_refinements = " << max_refinements <<endl;
    //Sets the numerical flux
    if (method=="upwind")
      update_flux=&upwind;
    else if(method != "lw")
    {
      cout <<"You incorrectly put method = "<<method<<endl;
      assert(false);
    }
    run_and_output(N_x, N_y, sigma_x, method, running_time,
                       initial_data_indicator, max_refinements);
}