//Solving Q_t + uQ_x + vQ_y = 0;

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
using namespace std;

//Returns true if real number is integer, false otherwise.
bool int_tester(double a)
{
    double dummy;
    double c = modf(a,&a);
    //cout << "Fractional part is "<<c<<endl;
    return (c<1e-4);
}

//Computes advection at x,y
void advection_velocity(double x, double y, double& u, double& v)
{
  u = 1.0;
  v = 1.0;
}

class Linear_Convection_2d
{
public:
    Linear_Convection_2d(const double N,
                         double lam_x, /*const double lam_y, */
                         string method, const double running_time,
                         int initial_data_indicator); 

    void run(bool output_indicator);
    void get_error(vector<double> &l1_vector, vector<double> &l2_vector, 
                   vector<double> &linfty_vector, vector<double> &snapshot_error);
private:
    void make_grid();
    void set_initial_solution();

    //Computes advection velocity at x,y
    void update_advection_velocity(double x, double y);
    //Computes flux_x(i+1/2,j), flux_y(i,j+1/2)
    void upwind_flux(int i, int j, double& flux_x,double& flux_y);
    void lw_flux(int i, int j, double& flux_x,double& flux_y);

    void upwind();
    void lw();

    double hat_function(double grid_point);
    double step_function(double grid_point); 
    double exp_func_25(double grid_point);
    double exp_func_100(double grid_point);


    void evaluate_error_and_output_solution(const int time_step_number,
                                            bool output_indicator);
    vector<double> grid_x,grid_y;

    double theta, u, v, xmin, xmax, ymin, ymax;

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

    double N, dx, dy, dt, t, running_time;
    double lam_x, lam_y, sigma_x, sigma_y;
    void update_ghost_values();
    void use_ghost_values();
    string method;
    int initial_data_indicator;
};

Linear_Convection_2d::Linear_Convection_2d(double N, 
                                           double lam_x,/* const double lam_y,*/
                                           string method,
                                           double running_time, 
                                           int initial_data_indicator):
                                           N(N), 
                                           lam_x(lam_x),
                                           running_time(running_time),
                                           method(method),
                                           initial_data_indicator(initial_data_indicator)
{
    theta = M_PI/4.0;
    u = 1.0, v = 1.0;
    //u = 1.0, v = 1.0;
    xmin = -1.0, xmax = 1.0, ymin = -1.0, ymax = 1.0;
    dx = (xmax - xmin) / (N), dy = (ymax-ymin)/(N);
    //In interval [0,1], if we run the loop for i = 0,1,...,n-1
    //and take the grid spacing to be 1/n, we won't reach the end of interval.
    t = 0.0;
    dt = lam_x*dx/(abs(u));
    /*//We want dt/u to be integer, this ensures it.
    dt = dt/u; Or maybe this step is not needed since dt already
    divided by u
    */
    lam_x = abs(u)*dt/dx;//lam_x updates
    sigma_x = u*dt/(dx), sigma_y = v*dt/(dy);
    //Since we have lam_x = |u|dt/dx,lam_y = |v|dt/dy, we can actually write
    //Thus, dt = lam_x*dx/|u|. Thus, lam_y = |v/u|*(dx/dy) * lam_x = 
    lam_y = abs(v/u)*dx/dy * lam_x;
    if (method == "lw"    &&
        (dx - dy < 1e-12) && 
        (dt/dx > 1.0/sqrt(u*u + v*v)))
    {
      cout << "WARNING - You are using lw with unstable sigma_x = "<< sigma_x <<endl;
    }
    snapshot_error = -1.0; //Put bad value for testing success.
    cout << "dx = " << dx << endl;
    cout << "dy = " << dy << endl;
    cout << "dt = " << dt << endl;
    cout << "lam_x = " <<lam_x << endl;
    cout << "lam_y = " <<lam_y << endl;
    //grid.resize(N);
    error.resize(N,N);
    grid_x.resize(N),grid_y.resize(N);
    initial_solution.resize(N,N);
    solution_old.resize(N,N,1);
    solution.resize(N,N,1);
    solution_exact.resize(N,N);
}

//Updates advection velocities at x,y, stored in class variables u,v
void Linear_Convection_2d::update_advection_velocity(double x, double y)
{
  u = 1.0, v = 1.0;
}

void Linear_Convection_2d::upwind_flux(int i, int j,
                                       double& flux_x, double& flux_y)
{
  //flux_x(i+1/2,j)
  flux_x = max(u,0.)*solution_old(i,j) + min(u,0.)*solution_old(i+1,j);
  //flux_y(i,j+1/2)
  flux_y = max(v,0.)*solution_old(i,j) + min(v,0.)*solution_old(i,j+1);
}

void Linear_Convection_2d::lw_flux(int i, int j, 
                                   double & flux_x, double& flux_y)
{
  
}

void Linear_Convection_2d::make_grid()
{
  //Note that you must run two for loops for a rectangular grid.
  for (unsigned int i = 0; i < N; i++)
  {
    grid_x[i] = (xmin+0.5*dx) + i * dx;
    grid_y[i] = (ymin+0.5*dy) + i * dy;
  }
}

void Linear_Convection_2d::update_ghost_values()
{
  //For readability, we do corners separately.
  solution_old(-1,-1)= solution_old(N-1,N-1);
  solution_old(N,N)  = solution_old(0,0);
  for (int j = 0; j<N;j++) //not doing corners
  {
    solution_old(-1,j) = solution_old(N-1,j);
    solution_old(N,j)  = solution_old(0,j);
  }
  //j = -1, N
  for (int i = 0; i<N;i++) 
  {
    solution_old(i,-1) = solution_old(i,N-1);
    solution_old(i,N)  = solution_old(i,0);
  }
  //We also set all ghost values in solution to equal zero.
  for (int k = -1; k<=N;k++)
  {
    solution(k,-1) = 0.0, solution(-1,k) = 0.0;
    solution(k,N)  = 0.0, solution(N,k) = 0.0;
  }
}

//This function adds anything that was added to ghost values to where it was
//supposed to go.
void Linear_Convection_2d::use_ghost_values()
{
  //Doing corners outside loop for readability
  solution(N-1,N-1) += solution(-1,-1), solution(0,0)+=solution(N,N);
  solution(-1,-1) = 0.0, solution(N,N) = 0;
  for (int j = 0; j<N;j++)
  {
    solution(N-1,j) += solution(-1,j), solution(0,j) += solution(N,j);
    solution(-1,j) = 0.0, solution(N,j) = 0.0;
  }
  for (int i = 0; i<N;i++)
  {
    solution(i,N-1) += solution(i,-1), solution(i,0)+=solution(i,N);
    solution(i,-1) = 0.0, solution(i,N)=0.0;
  }
}

void Linear_Convection_2d::set_initial_solution()
{
  double x,y;
  for (unsigned int i = 0; i < N; i++)
    for (unsigned int j = 0; j < N; j++)
    {
      x = (xmin+0.5*dx) + i*dx, y = (ymin+0.5*dy) + j*dy;
      switch (initial_data_indicator)
      {
      case 0:
        solution(i,j) = sin(2.0 * M_PI * x / (xmax - xmin))
                        * sin(2.0 * M_PI * y / (ymax - ymin));
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
      case 4:
        solution(i,j) = exp_func_100(x)*exp_func_100(y);
        break;
      default:
        cout << "You entered the wrong initial_data_indicator ";
        assert(false);
      }
      initial_solution(i,j) = solution(i,j); //Stored only
      //for snapshot error
    }
}

void Linear_Convection_2d::upwind()
{
  double x,y;
  double flux_x,flux_y; //flux_x(i+1/2,j), flux_y(i,j+1/2)
  //This loop computes the fluxes and adds them to where they are needed
  
  solution = 0.0;
  update_ghost_values();
  //We'd do solution = solution_old - dt/dx * (f_x(i+1/2,j)-f_x(i-1/2,j))
  //                                - dt/dx * (f_y(i,j+1/2)-f_y(i,j-1/2))
  for (int i = 0; i < N; i++) 
    for (int j = 0;j< N; j++)
    {
      x = (xmin + 0.5*dx) + i*dx, y = (ymin + 0.5*dy) + j*dy;
      update_advection_velocity(x,y); //Updates u,v
      upwind_flux(i,j,flux_x,flux_y); //flux_x(i+1/2,j), flux_y(i,j+1/2)
      //computed and stored in variables flux_x,flux_y.
      //Recall that, in FVM, solution updates as
      //Q_{i,j}^{n+1} = Q_{i,j}^n-(flux_x(i+1/2,j)-flux_x(i-1/2,j))*(dt/dx)
      //                         -(flux_y(i,j+1/2)-flux_y(i,j-1/2))*(dt/dy)
      //Q_{i+1,j}^{n+1} = Q_{i+1,j}^n-(flux_x(i+3/2,j)-flux_x(i+1/2,j))*(dt/dx)
      //                         -(flux_y(i+1,j+1/2)-flux_y(i+1,j-1/2))*(dt/dy)
      //Q_{i,j+1}^{n+1} = Q_{i,j+1}^n-(flux_x(i+1/2,j+1)-flux_x(i-1/2,j+1))*(dt/dx)
      //                         -(flux_y(i,j+3/2)-flux_y(i,j+1/2))*(dt/dy)
      //Q_{i+1,j+1}^{n+1} = Q_{i+1,j+1}^n-(flux_x(i+3/2,j+1)-flux_x(i+1/2,j+1))*(dt/dx)
      //                         -(flux_y(i+1,j+3/2) - flux_y(i+1,j+1/2))*(dt/dx)
      //Q_{i-1,j-1}^{n+1} = Q_{i-1,j-1}^n-(flux_x(i-1/2,j-1)-flux_x(i-3/2,j-1))*(dt/dx)
      //                                  -(flux_y(i-1,j-1/2)-flux_y(i-1,j-3/2))*(dt/dy)
      solution(i,j)     += -flux_x*(dt/dx) - flux_y*(dt/dy);
      solution(i,j+1)   +=  flux_y*(dt/dy);
      solution(i+1,j)   +=  flux_x*(dt/dx);
      solution(i,j)     +=  solution_old(i,j);  
    }
    use_ghost_values();
}

void Linear_Convection_2d::lw()
{
  double x,y;
  double flux_x,flux_y; //flux_x(i+1/2,j), flux_y(i,j+1/2)
  //This loop computes the fluxes and adds them to where they are needed
  
  solution = 0.0;
  update_ghost_values();
  //We'd do solution = solution_old - dt/dx * (f_x(i+1/2,j)-f_x(i-1/2,j))
  //                                - dt/dx * (f_y(i,j+1/2)-f_y(i,j-1/2))
  for (int i = 0; i < N; i++) 
    for (int j = 0;j< N; j++)
    {
      x = (xmin + 0.5*dx) + i*dx, y = (ymin + 0.5*dy) + j*dy;
      update_advection_velocity(x,y); //Updates u,v
      lw_flux(i,j,flux_x,flux_y); //flux_x(i+1/2,j), flux_y(i,j+1/2)
      solution(i,j+1)   +=  flux_y*(dt/dy);
      solution(i+1,j)   +=  flux_x*(dt/dx);
      solution(i,j)     +=  solution_old(i,j);  
    }
    use_ghost_values();
}

void Linear_Convection_2d::evaluate_error_and_output_solution(int time_step_number,bool output_indicator)
{
  double x,y;
  for (unsigned int i = 0; i < N; i++)
    for (unsigned int j = 0; j < N; j++)
    {
      x = (xmin+0.5*dx) + i*dx, y = (ymin+0.5*dy) + j*dy;
      switch (initial_data_indicator)
      {
      case 0:
        solution_exact(i,j) = sin(2.0 * M_PI * (x - u*t) / (xmax - xmin))
                              * sin(2.0 * M_PI * (y - v*t) / (ymax - ymin));
        break;
      case 1:
        solution_exact(i,j) = hat_function(x - u*t)
                              * hat_function(y - v*t);
        break;
      case 2:
        solution_exact(i,j) = step_function(x - u*t)
                              * step_function(y - v*t);
        break;
      case 3:
        solution_exact(i,j) = exp_func_25(x-u*t)
                              * exp_func_25(y-v*t);
        break;
      case 4:
        solution_exact(i,j) = exp_func_100(x-u*t)
                              * exp_func_100(y-v*t);
        break;
      default:
          cout << "You entered the wrong initial_data_indicator ";
          assert(false);
      }
    }
  if (output_indicator==true && time_step_number%5==0)
  {
  vtk_anim_sol(grid_x,grid_y,
        solution,
        t, time_step_number/5,
        "approximate_solution");
    vtk_anim_sol(grid_x,grid_y,
        solution_exact,
        t, time_step_number/5,
        "exact_solution");
  }
  for (unsigned int i = 0; i < N; i++)
    for (unsigned int j = 0; j < N; j++)
    {
      error(i,j) = abs(solution(i,j) - solution_exact(i,j));
    }
}

void Linear_Convection_2d::run(bool output_indicator)
{
  make_grid();
  int time_step_number = 0; 
  set_initial_solution(); //sets solution to be the initial data
  evaluate_error_and_output_solution(time_step_number,output_indicator);
  while (t < running_time) //compute solution at next time step using solution_old
  {        
    solution_old = solution;//update solution_old for next time_step
    //This requires puttingg ghost cells in solution
    //We should update array2d.h to change this. 
    if (method == "upwind")
      upwind();
    else
      assert(false);
    time_step_number += 1;
    t = t + dt; 
    evaluate_error_and_output_solution(time_step_number, output_indicator);
  }
  cout << "For N = " << N<<", we took ";
  cout << time_step_number << " steps." << endl;
}

void run_and_output(double N, double cfl,
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
    Linear_Convection_2d solver(N, cfl, method, running_time,
                                initial_data_indicator);
    //We calculate time takenṣ in our refinement.
    struct timeval begin, end; 
    gettimeofday(&begin, 0);
    solver.run(refinement_level==max_refinements-1);//Output only last soln
    gettimeofday(&end, 0);
    long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
    double elapsed = seconds + microseconds * 1e-6;
    cout << "Time taken by this iteration is " << elapsed << " seconds." << endl;
    solver.get_error(l1_vector,l2_vector,linfty_vector,snapshot_vector);//push_back resp. error.
    error_vs_h << 2.0/N << " " << linfty_vector[refinement_level] << "\n";
    N = 2.0 * N;
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

double Linear_Convection_2d::interval_part(double x)
{
  if (x > xmax)
    return x - ceil((x - xmax) / (xmax - xmin)) * (xmax - xmin);
  else if (x < xmin)
    return x + ceil((xmin - x) / (xmax - xmin)) * (xmax - xmin);
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
  if (grid_point < xmin + (xmax - xmin) / 4.0 || 
    grid_point > xmax - (xmax - xmin) / 4.0)
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

double Linear_Convection_2d::exp_func_100(double grid_point)
{
  grid_point = interval_part(grid_point);
  return exp(-100.0*(grid_point-0.5)*(grid_point-0.5));
}

void Linear_Convection_2d::get_error(vector<double> &l1_vector,
                                     vector<double> &l2_vector,
                                     vector<double> &linfty_vector,
                                     vector<double> &snapshot_vector)
{
    double l1 = 0.0;
    double l2 = 0.0;
    double linfty = 0.0;
    for (int j = 0; j < N; j++) 
      for (int i = 0; i < N; i++)
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
    if (int_tester(0.5*running_time*u) == true &&
        int_tester(0.5*running_time*v) == true)
    {
      cout << "Snapshot error successfuly computed.\n";
      for (unsigned int i = 0; i < N; i++)
        for (unsigned int j = 0; j < N; j++)
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
    if (argc != 6)
    {
      cout << "Incorrect format, use" << endl;
      cout << "./fd2d upwind/lw/ct_upwind(for the respective method)";
      cout << " sigma_x running_time initial_data_indicator " ;
      cout << "max_refinements" << endl;
      cout << "Choices for method"<<endl;
      cout << "upwind, lw, ct_upwind, m_roe"<<endl;
      cout << "Choices for initial data "<<endl;
      cout << "0 - smooth_sine"<<endl<<"1 - hat"<<endl;
      cout << "2 - step"<<endl;
      cout << "3 - exp_func_25" <<endl;
      assert(false);
    }
    string method = argv[1];
    cout << "method = " << method << endl;
    double N = 30.0;
    double sigma_x = stod(argv[2]);
    cout << "sigma_x = " << sigma_x << endl;
    double running_time = stod(argv[3]);
    cout << "running_time = " << running_time << endl;
    int initial_data_indicator = stoi(argv[4]);
    cout << "initial_data_indicator = " << initial_data_indicator << endl;
    unsigned int max_refinements = stoi(argv[5]);
    cout << "max_refinements = " << max_refinements <<endl;
    run_and_output(N, sigma_x, method, running_time,
                       initial_data_indicator, max_refinements);
}
