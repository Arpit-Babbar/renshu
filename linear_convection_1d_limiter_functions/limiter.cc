#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <string>
#include <stdio.h>
#include <sys/time.h>

using namespace std;

//Defining y = x1 + a2*x2
void add(const vector<double> &x1,
         const double a2, const vector<double> &x2,
         vector<double> &y)
{
    if (x1.size()!=x2.size() || x2.size() != y.size())
        {
            cout << "Tried y = x1 + a2*x2 with inappropriately sized vectors.";
            cout << "Size of chosen x1 = "<< x1.size() <<". ";
            cout << "Size of chosen x2 = "<< x2.size()<<". ";
            cout << "Size of chosen y = "<< y.size()  <<". ";
            assert(false);
        }
    for (unsigned int i = 0; i < y.size(); i++)
    {
        y[i] = x1[i] + a2*x2[i];
    }
}

//Defining y = x1 + a2*x2 + a3*x3
void add(const vector<double> &x1,
         const double a2, const vector<double> &x2,
         const double a3, const vector<double> &x3,
         vector<double> &y)
{
    if (x1.size()!=x2.size() || x2.size() != x3.size() || x3.size() != y.size())
        {
            cout << "Tried y = x1 + a2*x2 + a3*x3 with inappropriately sized vectors.";
            cout << "Size of chosen x1 = "<< x1.size() <<". ";
            cout << "Size of chosen x2 = "<< x2.size() <<". ";
            cout << "Size of chosen x3 = "<< x3.size() <<". ";
            cout << "Size of chosen y = "<< y.size()   <<". ";
            assert(false);
        }
    for (unsigned int i = 0; i < y.size(); i++)
    {
        y[i] = x1[i] + a2*x2[i] + a3*x3[i];
    }
}

void output_vectors_to_file(string file_name, vector<double> &grid,
                            vector<double> &solution_old, vector<double> &solution_exact)
{
    ofstream output_solution;
    output_solution.open(file_name);
    for (unsigned int j = 0; j < grid.size(); j++)
    {
        output_solution << grid[j] << " " << solution_old[j] << " " << solution_exact[j] << "\n";
    }
    output_solution.close();
}

void output_vectors_to_file(string file_name, vector<double> &grid, vector<double> &solution_old)
{
    ofstream output_solution;
    output_solution.open(file_name);
    for (unsigned int j = 0; j < grid.size(); j++)
    {
        output_solution << grid[j] << " " << solution_old[j] << "\n";
    }
    output_solution.close();
}

class Linear_Convection_1d
{
public:
    Linear_Convection_1d(const double n_points, const double cfl,
                         string method, const double running_time,
                         int initial_data_indicator); //input parameters
    //and which method to use - Lax-Wendroff or RK4

    void run(); //true when output is to be given, and false when it doesn't;
    void output_final_error();
    void get_error(vector<double> &l1_vector, vector<double> &l2_vector, vector<double> &linfty_vector);
private:
    void make_grid(); // Grid is needed for writing output, and defining exact solution.
    void set_initial_solution();
    //removed initial_data vector, as it was wasteful
    //We'd use this to compute the rk4 slopes k_i's and temporarily store them in solution_new.
    //More precisely, it does solution = solution_old + factor * u
    void lax_wendroff();  
    void foup(); //First order upwind scheme
    
    void rhs_function();
    
    void rk2_solver();
    void rk3_solver();

    void rhs_soup();
    void rhs_soup_minmod(); //This gives RHS of the system of ODEs and stores it to where rhs points.
    

    void soup_rk2(); //Second order upwind Scheme
    void soup_rk3();
    void soup_rk2_minmod(); //Second order upwind Scheme with min-mod limiter
    void soup_rk3_minmod(); //Second order upwind Scheme with min-mod limiter

    double minmod_limiter(int i); //psi(r_i)

    double hat_function(double grid_point);
    double step_function(double grid_point); //Functions for initial data.
                                             //on which we apply RK4, and stores it in k.
    double sine_wave(double grid_point);

    void evaluate_error_and_output_solution(const int time_step_number);

    double coefficient, x_min, x_max;

    double interval_part(double); //If a point is not in the interval [x_min,x_max], 
    //this function sends it there by translating in x_max-x_min sized steps.

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

    string method;
    int initial_data_indicator;
};

Linear_Convection_1d::Linear_Convection_1d(double n_points, double cfl, string method,
                                           double running_time, int initial_data_indicator):
                                           n_points(n_points), cfl(cfl),
                                           running_time(running_time), method(method),
                                           initial_data_indicator(initial_data_indicator)
{
    coefficient = 1.0;
    sigma = cfl * coefficient / abs(coefficient); //sigma = coefficient*dt/h
    x_min = -1.0, x_max = 1.0;
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

void Linear_Convection_1d::make_grid()
{
    for (unsigned int i = 0; i < n_points; i++)
    {
        grid[i] = x_min + i * h;
    }
}

void Linear_Convection_1d::set_initial_solution()
{
    for (unsigned int i = 0; i < n_points; i++)
    {
        switch (initial_data_indicator)
        {
        case 0:
            solution[i] = sin(2.0 * M_PI * grid[i] / (x_max - x_min));
            break;
        case 1:
            solution[i] = hat_function(grid[i]);
            break;
        case 2:
            solution[i] = step_function(grid[i]);
            break;
        case 3:
            solution[i] = sine_wave(grid[i]); 
            //sine function supported in [-0.5,0.5] zero in rest of [-1,1].
            break;
        default:
            cout << "You entered the wrong initial_data_indicator ";
            assert(false);
        }
    }
}

//This computes the rhs of the system of ODEs on which we apply rk4.

void Linear_Convection_1d::lax_wendroff()
{
    solution[0] = solution_old[0] 
                     - 0.5 * cfl * (solution_old[1] - solution_old[n_points - 1])
                     + 0.5 * cfl * cfl * (solution_old[n_points - 1] 
                                          - 2.0 * solution_old[0] + solution_old[1]);
    //For easy readability, whenever there is a line break in an ongoing bracket,
    //the next line starts from where the bracket opens.
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

void Linear_Convection_1d::foup() //first order upwind scheme
{
  solution[0] = (1.0 - sigma) * solution_old[0] + cfl * solution_old[n_points - 1];
  for (int i = 1; i < n_points; i++)
  {
    solution[i] = (1.0 - sigma) * solution_old[i] + cfl * solution_old[i - 1];
  }
}

double Linear_Convection_1d::minmod_limiter(int i) //Psi(r_i)
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
    return min(r,1.0);
    }
}


void Linear_Convection_1d::rhs_function()
{
  if (method == "soup_rk2" || method == "soup_rk3")
  {
    (*rhs)[0] = - (coefficient/h) * (1.5 * solution_old[0] 
                                 - 2.0 * solution_old[n_points-1] 
                                 + 0.5 * solution_old[n_points-2]);
    (*rhs)[1] = - (coefficient/h) * (1.5 * solution_old[1] 
                                 - 2.0 * solution_old[0] 
                                 + 0.5 * solution_old[n_points-1]);
    for (int i = 2; i < n_points; i++)
    {
        (*rhs)[i] = - (coefficient/h) * (1.5 * solution_old[i] 
                                    - 2.0 * solution_old[i-1] 
                                    + 0.5 * solution_old[i-2]);
    }
  }
  else if (method == "soup_rk3_minmod" || method == "soup_rk2_minmod")
  {
    (*rhs)[0] = -(coefficient/h) * ( (solution_old[0] - solution_old[n_points-1])
                                  +0.5 * minmod_limiter(0)
                                        * (solution_old[0] - solution_old[n_points-1])
                                  -0.5 * minmod_limiter(n_points-1)
                                        * (solution_old[n_points-1] - solution_old[n_points-2]));
    (*rhs)[1] = -(coefficient/h) * ( (solution_old[1] - solution_old[0])
                                  +0.5 * minmod_limiter(1)
                                        * (solution_old[1] - solution_old[0])
                                  -0.5 * minmod_limiter(0)
                                        * (solution_old[0] - solution_old[n_points-1]));
    for (int i = 2; i<n_points; i++)
    {
      (*rhs)[i] = -(coefficient/h) * ( (solution_old[i] - solution_old[i-1])
                                  +0.5 * minmod_limiter(i)
                                        * (solution_old[i] - solution_old[i-1])
                                  -0.5 * minmod_limiter(i-1)
                                        * (solution_old[i-1] - solution_old[i-2]));
    }
  }
  else 
  {
    cout << "Incorrect method in rhs_function(method)" << endl;
    assert(false);
  }
}


void Linear_Convection_1d::rk3_solver()
{
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
}

void Linear_Convection_1d::rk2_solver()
{
    rhs = &temp;
    //k1 = rhs_function(solution_old)
    rhs_function(); //depends on method
    //Temporarily putting u^{n+1} = solution = solution_old + 0.5 * dt * k1
    add(solution_old, 0.5 * dt, temp, solution);
    //Next, putting k1 = rhs_function(solution) = rhs_function(solution+old + 0.5)
    rhs_function(); 
    add(solution_old,dt, temp, solution);
}
void Linear_Convection_1d::evaluate_error_and_output_solution(int time_step_number)
{
    for (unsigned int i = 0; i < n_points; i++)
    {
        switch (initial_data_indicator)
        {
        case 0:
            solution_exact[i] = sin(2.0 * M_PI * (grid[i] - coefficient * t) / (x_max - x_min));
            break;
        case 1:
            solution_exact[i] = hat_function(grid[i] - coefficient * t);
            break;
        case 2:
            solution_exact[i] = step_function(grid[i] - coefficient * t);
            break;
        case 3:
            solution_exact[i] = sine_wave(grid[i] - coefficient * t);
            break;
        default:
            cout << "Wrong initial_data_indicator entered.";
            assert(false);
        }
    }
    
    for (unsigned int j = 0; j < n_points; j++)
    {
        error[j] = abs(solution[j] - solution_exact[j]);
    }
    string solution_file_name = "solution_";
    solution_file_name += to_string(time_step_number) + ".txt";
    output_vectors_to_file(solution_file_name, grid, solution, solution_exact);
}

void Linear_Convection_1d::output_final_error()
{
    string error_file_name = "finalerror_";
    error_file_name += method + ".txt";
    output_vectors_to_file(error_file_name, grid, error);
}

void Linear_Convection_1d::run()
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
        else if (method == "soup_rk3" ||  method == "soup_rk3_minmod")
            rk3_solver();
            //rk3_solver in which the rhs would depend on method
        else if (method == "soup_rk2" || method == "soup_rk2_minmod")
            rk2_solver();
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

    for (unsigned int refinement_level = 0; refinement_level <= max_refinements;
        refinement_level++)
    {
        Linear_Convection_1d solver(n_points, cfl, method, running_time,
                                    initial_data_indicator);
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

double Linear_Convection_1d::interval_part(double x)
{
    if (x > x_max)
        return x - ceil((x - x_max) / (x_max - x_min)) * (x_max - x_min);
    else if (x < x_min)
        return x + ceil((x_min - x) / (x_max - x_min)) * (x_max - x_min);
    else
        return x;
}

double Linear_Convection_1d::hat_function(double grid_point)
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

double Linear_Convection_1d::step_function(double grid_point)
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

double Linear_Convection_1d::sine_wave(double grid_point)
{
    double value;
    grid_point = interval_part(grid_point);
    if (grid_point < x_min + (x_max - x_min) / 4.0 ||
        grid_point > x_max - (x_max - x_min) / 4.0)
        value = 0.0;
    else
        value = sin(4 * M_PI * grid_point / (x_max - x_min));
    return value;
}

void Linear_Convection_1d::get_error(vector<double> &l1_vector, 
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
    double n_points = 50.0;
    double cfl = stod(argv[2]);
    cout << "cfl = " << cfl << endl;
    double running_time = stod(argv[3]);
    cout << "running_time = " << running_time << endl;
    int initial_data_indicator = stoi(argv[4]);
    cout << "initial_data_indicator = " << initial_data_indicator << endl;
    unsigned int max_refinements = stoi(argv[5]);
    cout << "max_refinements = " << max_refinements <<endl;
    run_and_get_output(n_points, cfl, method, running_time, initial_data_indicator,
                       max_refinements);
}