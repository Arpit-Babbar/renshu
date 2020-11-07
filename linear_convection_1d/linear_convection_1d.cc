#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <string>
#include <stdio.h>
#include <sys/time.h>

using namespace std;


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
    void temporary_update_solution(const double factor, const vector<double> &u);
    //We'd use this to compute the rk4 slopes k_i's and temporarily store them in solution_new.
    //More precisely, it does k = solution_old + factor * u
    void temporary_update_solution(const double factor0, const vector<double> &u0,
                                   const double factor1, const vector<double> &u1);
    void rk4_solver(); //Gives the solution at next time step using RK4
    void rk3_solver();
    void rk2_solver();
    void lax_wendroff();                                           
    void rhs_function(const vector<double> &u, vector<double> &k); //This gives RHS of the system of ODEs
    double hat_function(double grid_point);
    double step_function(double grid_point); //Functions for initial data.
                                             //on which we apply RK4, and stores it in k.

    void evaluate_error_and_output_solution(const int time_step_number);

    double coefficient, x_min, x_max;

    double interval_part(double); //If a point is not in the interval [x_min,x_max], 
    //this function sends it there by translating in x_max-x_min sized steps.

    vector<double> grid;

    vector<double> solution_old; //Solution at previous step
    vector<double> solution; //Solution at present step

    vector<double> solution_exact; //Exact solution at present time step

    vector<double> error; //This will store the maximum error at a grid point in all time-steps
    double l2_error;

    //h denotes the spatial distance between grid points
    double n_points, h, dt,t, cfl, running_time; //h = 1/n_points just included for easy typing

    //slopes needed by rk4
    vector<double> k1, k2, k3, k4;

    //This computes solution_new = u^{n+1} = u^n + dt/6 * (k1 + 2.0*k2 + 2.0*k3 + k4)
    void compute_solution_new_using_ki(int);

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
    k1.resize(n_points);
    k2.resize(n_points);
    k3.resize(n_points);
    k4.resize(n_points);
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
            solution_old[i] = sin(2.0 * M_PI * grid[i] / (x_max - x_min));
            break;
        case 1:
            solution_old[i] = hat_function(grid[i]);
            break;
        case 2:
            solution_old[i] = step_function(grid[i]);
            break;
        default:
            cout << "You entered the wrong initial_data_indicator ";
            assert(false);
        }
    }
}

//This computes the rhs of the system of ODEs on which we apply rk4.
void Linear_Convection_1d::rhs_function(const vector<double> &u, vector<double> &k)
{
    if (u.size() != k.size() || k.size() != solution_old.size())
    {
        cout << "You have used rhs_function to do rhs_funciton(&u,&k)"
             << " with inappropriately sized u,k. Size of your u is " << u.size();
        cout << ", size of k is " << k.size()
             << ". Both the sizes should equal " << solution_old.size() << endl;
        assert(false);
    }
    k[0] = -(u[1] - u[n_points - 1]) / (2.0 * h); //left end point
    for (int j = 1; j < n_points - 1; j++)
    {
        k[j] = -(u[j + 1] - u[j - 1]) / (2.0 * h);
    }
    k[n_points - 1] = -(u[0] - u[n_points - 2]) / (2.0 * h); //right end point
}

//k = solution_old + factor*u
void Linear_Convection_1d::temporary_update_solution(const double factor, const vector<double> &k0)
{
    if (k0.size() != solution.size() || solution.size() != solution_old.size())
    {
        cout << "You have used temporary_update_solution to do 'k = solution_old + factor * k0 ";
        cout << "with inappropriately sized k0,k. Size of your k0 is " << k0.size();
        cout << ", size of k is " << solution.size()
             << ". Both the sizes should equal " << solution_old.size() << endl;
        assert(false);
    }
    for (unsigned int i = 0; i < solution.size(); i++)
    {
        solution[i] = solution_old[i] + factor * k0[i];
    } //k = solution_old + factor * k0
}
void Linear_Convection_1d::temporary_update_solution(const double factor0, const vector<double> &k0,
                                                     const double factor1, const vector<double> &k1)
{
    if (k0.size() != solution.size() || solution.size() != solution_old.size())
    {
        cout << "You have used temporary_update_solution to do 'k = solution_old + factor * k0 ";
        cout << "with inappropriately sized k0,k. Size of your k0 is " << k0.size();
        cout << ", size of k is " << solution.size();
        cout << ". Both the sizes should equal " << solution_old.size() << endl;
        assert(false);
    }
    for (unsigned int i = 0; i < solution.size(); i++)
    {
        solution[i] = solution_old[i] + factor0 * k0[i] + factor1 * k1[i];
    } //k = solution_old + factor * k0
}

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

void Linear_Convection_1d::rk4_solver()
{
    //This innocent step was what was causing trouble
    /* solution_old = solution; */

    //k1 = rhs_function(solution_old)
    rhs_function(solution_old, k1);
    //To save memory, we are going to temporarily store some things in solution_new
    //, which is supposed to be u^{n+1}
    //Temporarily putting 'u^{n+1}' = solution_new = solution_old + dt/2 * k1
    temporary_update_solution(dt / 2.0, k1);
    //So, computing k2 = rhs_function(u^n + dt/2 * k1)
    rhs_function(solution, k2);
    //Similarly, temporarily putting 'u^{n+1}'=solution_new = u^n + dt/2 *k2
    temporary_update_solution(dt / 2.0, k2);
    //Computing k3 = rhs_function(solution_old + dt/2 *k2)
    rhs_function(solution, k3);
    //Temporarily putting 'u^{n+1}' = solution_new = u^n + dt * k3
    temporary_update_solution(dt, k3);
    //Computing k4 = rhs_function(solution_old + dt *k3)
    rhs_function(solution, k4);

    //This computes solution_new = u^{n+1} = u^n + dt/6 * (k1 + 2.0*k2 + 2.0*k3 + k4)
    for (unsigned int i = 0; i < n_points; i++)
        solution[i] = solution_old[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}

void Linear_Convection_1d::rk3_solver()
{
    //k1 = rhs_function(solution_old)
    rhs_function(solution_old, k1);
    //Temporarily putting u^{n+1} = solution_new = solution_old + dt/2 * k1
    temporary_update_solution(dt / 2.0, k1);
    //So, computing k2 = rhs_function(u^n + dt/2 * k1)
    rhs_function(solution, k2);
    //Similarly, temporarily putting u^{n+1}=solution_new = u^n -k1 * dt + 2 k2 * dt
    temporary_update_solution(-dt, k1, 2.0 * dt, k2);
    //Computing k3 = rhs_function(solution_old - dt k1 + 2 dt k2)
    rhs_function(solution, k3);

    //This computes solution_new = u^{n+1} = u^n + dt/6 * (k1 + 2.0*k2 + 2.0*k3 + k4)
    for (unsigned int i = 0; i < n_points; i++)
        solution[i] = solution_old[i] + dt / 6.0 * (k1[i] + 4.0 * k2[i] + k3[i]);
}

void Linear_Convection_1d::rk2_solver()
{
    //k1 = rhs_function(solution_old)
    rhs_function(solution_old, k1);
    //Temporarily putting u^{n+1} = solution_new = solution_old + dt * k1
    temporary_update_solution(dt, k1);
    //So, computing k2 = rhs_function(u^n + dt * k1)
    rhs_function(solution, k2);

    //This computes solution_new = u^{n+1} = u^n + dt/6 * (k1 + 2.0*k2 + 2.0*k3 + k4)
    for (unsigned int i = 0; i < n_points; i++)
        solution[i] = solution_old[i] + dt / 2.0 * (k1[i] + k2[i]);
}

void Linear_Convection_1d::evaluate_error_and_output_solution(int time_step_number)
{
    for (unsigned int j = 0; j < n_points; j++)
    {
        switch (initial_data_indicator)
        {
        case 0:
            solution_exact[j] = sin(2.0 * M_PI * (grid[j] - t) / (x_max - x_min));
            break;
        case 1:
            solution_exact[j] = hat_function(grid[j] - t);
            break;
        case 2:
            solution_exact[j] = step_function(grid[j] - t);
            break;
        default:
            cout << "You entered the wrong initial_data_indicator ";
            assert(false);
        }
    }
    
    for (unsigned int j = 0; j < n_points; j++)
    {
        error[j] = abs(solution[j] - solution_exact[j]);
    }
    string solution_file_name = "solution_";
    solution_file_name += to_string(time_step_number) + "_" + method + ".txt";
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
    set_initial_solution(); //sets solution_old to be the initial data
    int time_step_number = 0; //The bug I had made here was that I was running the solver
    //even for time_step_number = 0, which is wrong.
    evaluate_error_and_output_solution(time_step_number);
    while (t < running_time) //compute solution at next time step using solution_old
    {
        if (method == "rk4")
            rk4_solver();
        else if (method == "rk3")
            rk3_solver();
        else if (method == "rk2")
            rk2_solver();
        else if (method == "lw")
            lax_wendroff();
        else
            assert(false);
        solution_old = solution;//update solution_old to be used at next time step
        time_step_number += 1;
        t = t + dt; //Need to update the time_step_number because the updated version
        //is needed by the next function. We could have done it outside of the while
        //loop, but that was causing some duplication. Still not sure what the 
        //optimum thing to do is.
        evaluate_error_and_output_solution(time_step_number);
    }
    cout << "For n_points = " << n_points<<", we took ";
    cout << time_step_number << " steps." << endl;
}

void run_and_get_output(double n_points, double cfl,
                        string method, double running_time,
                        int initial_data_indicator,
                        int max_refinements)
{
    ofstream error_vs_h;
    error_vs_h.open("error_vs_h.txt");
    Linear_Convection_1d solver(n_points, cfl, method, running_time, initial_data_indicator);
    solver.run();
    vector<double> linfty_vector;
    vector<double> l2_vector;
    vector<double> l1_vector;
    int refinement_level = 0;

    while (refinement_level <= max_refinements)
    {
        solver.get_error(l1_vector,l2_vector,linfty_vector);//push_back resp. error.
        error_vs_h << 1.0/n_points << " " << linfty_vector[refinement_level] << "\n";
        n_points = 2.0 * n_points;
        solver = Linear_Convection_1d(n_points, cfl, method,
                                      running_time, initial_data_indicator);
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
        struct timeval begin, end;
        gettimeofday(&begin, 0);
        solver.run();
        gettimeofday(&end, 0);
        long seconds = end.tv_sec - begin.tv_sec;
        long microseconds = end.tv_usec - begin.tv_usec;
        double elapsed = seconds + microseconds * 1e-6;
        cout << "Time taken by this iteration is " << elapsed << " seconds." << endl;
        refinement_level++;
    }
    error_vs_h.close();
    solver.output_final_error();
    cout << "After " << refinement_level << " refinements, l_infty error = ";
    cout << linfty_vector[max_refinements-1] << endl;
    cout << "The L2 error is " << l2_vector[max_refinements-1] << endl;
    cout << "The L1 error is " << l1_vector[max_refinements-1] << endl;
    cout << endl;         
}

double Linear_Convection_1d::interval_part(double x)
{
    if (x > x_max)
        return x - ceil((x_max - x_min) * (x - x_max)) / (x_max - x_min);
    else if (x < x_min)
        return x + ceil((x_max - x_min) * (x_min - x)) / (x_max - x_min);
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
    if (grid_point < x_min + (x_max - x_min) / 4.0 || grid_point > x_max - (x_max - x_min) / 4.0)
        value = 0.0;
    else
        value = 1.0;
    return value;
}
void Linear_Convection_1d::get_error(vector<double> &l1_vector, vector<double> &l2_vector,
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
        cout << "./output lw/rk4/rk3/rk2(for the respective method) cfl running_time 0/1/2(for " ;
        cout << "sin/hat/discts initial data) max_refinements" << endl;
        assert(false);
    }
    string method = argv[1];
    cout << "method = " << method << endl;
    double n_points = 15.0;
    double cfl = stod(argv[2]);
    cout << "cfl = " << cfl << endl;
    double running_time = stod(argv[3]);
    cout << "running_time = " << running_time << endl;
    int initial_data_indicator = stoi(argv[4]);
    cout << "initial_data_indicator = " << initial_data_indicator << endl;
    int max_refinements = stoi(argv[5]);
    cout << "max_refinements = " << max_refinements <<endl;
    run_and_get_output(n_points, cfl, method, running_time, initial_data_indicator, max_refinements);
}
