#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include<string>

#include <algorithm>
#include <functional> //Used to define addition of vectors

using namespace std;
class Linear_Convection_1d
{
public:
    Linear_Convection_1d(const double n_points, const double cfl, string method, const double running_time, int initial_data_indicator); //input parameters
    //and which method to use - Lax-Wendroff or RK4
    void run_and_output_results();

private:
    void make_grid(); // Grid is needed for writing output, and defining exact solution.
    void set_initial_data();
    void temporary_update_solution(const double factor, const vector<double> &u, vector<double> &k);
    //We'd use this to compute the rk4 slopes k_i's and temporarily store them in solution_new.
    //More precisely, it does k = solution_old + factor * u
    void rk4_solver();                                             //Gives the solution at next time step using RK4
    void lax_wendroff();                                           //Gives the solutions at next time step using Lax-Wendroff
    void rhs_function(const vector<double> &u, vector<double> &k); //This gives RHS of the system of ODEs
    double hat_function(double grid_point);
    double step_function(double grid_point); //Functions for initial data.
    //on which we apply RK4, and stores it in k.

    void output_error_and_results(const int time_step_number);

    double coefficient = 1.0;

    double x_min = - 1.0;
    double x_max = 1.0;

    vector<double> grid;

    vector<double> solution_old; //Solution at previous step
    vector<double> solution_new; //Solution at present step

    vector<double> solution_exact; //Exact solution at present time step

    vector<double> initial_data;

    vector<double> error;

    //h denotes the spatial distance between grid points
    double n_points, h, dt, cfl, running_time; //h = 1/n_points just included for easy typing

    //slopes needed by rk4
    vector<double> k1, k2, k3, k4;

    //This computes solution_new = u^{n+1} = u^n + dt/6 * (k1 + 2.0*k2 + 2.0*k3 + k4)
    void compute_solution_new_using_ki();

    string method;
    int initial_data_indicator;
};

Linear_Convection_1d::Linear_Convection_1d(double n_points, double cfl, string method, double running_time, int initial_data_indicator) : n_points(n_points), cfl(cfl), running_time(running_time), method(method), initial_data_indicator(initial_data_indicator)
{
    h = (x_max - x_min) / n_points;
    dt = cfl * h / abs(coefficient);
    cout << "Spatial grid points are at gap h = " << h << endl;
    cout << "Time step dt = " << dt << endl;
    grid.resize(n_points);
    initial_data.resize(n_points);
    error.resize(n_points);
    solution_old.resize(n_points);
    solution_new.resize(n_points);
    solution_exact.resize(n_points);
    cout << "Number of spatial grid points is " << n_points << endl;
    k1.resize(n_points);
    k2.resize(n_points);
    k3.resize(n_points);
    k4.resize(n_points);
};

void Linear_Convection_1d::make_grid()
{
    for (int i = 0; i < n_points; i++)
    {
        grid[i] = x_min + i * h;
    };
};

double Linear_Convection_1d::hat_function(double grid_point)
{
    double value;
    if (grid_point < -1.0/2.0 || grid_point > 1.0/2.0)
        value = 0;
    else if (grid_point <= 0.0)
        value = grid_point + 1.0/2.0;
    else if (grid_point >= 0.0)
        value = 1.0/2.0 - grid_point;
    return value;
}

void Linear_Convection_1d::set_initial_data()
{
    switch (initial_data_indicator)
    {
    case 0:
        for (int i = 0; i < n_points; i++)
        {
            initial_data[i] = sin(/*2 * M_PI **/ grid[i]);
        }
        break;
    case 1:
        for (int i = 0; i < n_points; i++)
        {
            initial_data[i] = hat_function(grid[i]);
        }
        break;
    case 2:
        for (int i = 0; i < n_points; i++)
        {
            initial_data[i] = step_function(grid[i]);
        }
        break;
    default:
        cout << "You entered the wrong initial_data_indicator ";
        assert(false);
    };
}

//This computes the rhs of the system of ODEs on which we apply rk4.
void Linear_Convection_1d::rhs_function(const vector<double> &u, vector<double> &k)
{
    if (u.size() != k.size() || k.size() != solution_old.size())
    {
        cout << "You have used rhs_function to do rhs_funciton(&u,&k)"
             << " with inappropriately sized u,k. Size of your u is " << u.size() << ", size of k is " << k.size()
             << ". Both the sizes should equal " << solution_old.size() << endl;
        assert(false);
    }
    k[0] = -coefficient * (u[1] - u[n_points - 1]) / (2 * h); //left end point
    for (int j = 1; j < n_points - 1; j++)
    {
        k[j] = -coefficient * (u[j + 1] - u[j - 1]) / (2 * h);
    }
    k[n_points - 1] = -coefficient * (u[0] - u[n_points - 2]) / (2 * h); //right end point
};

//k = solution_old + factor*u
void Linear_Convection_1d::temporary_update_solution(const double factor, const vector<double> &k0, vector<double> &k)
{
    if (k0.size() != k.size() || k.size() != solution_old.size())
    {
        cout << "You have used temporary_update_solution to do 'k = solution_old + factor * k0 "
             << "with inappropriately sized k0,k. Size of your k0 is " << k0.size() << ", size of k is " << k.size()
             << ". Both the sizes should equal " << solution_old.size() << endl;
        assert(false);
    }
    for (int i = 0; i < k.size(); i++)
    {
        k[i] = solution_old[i] + factor * k0[i];
    } //k = solution_old + factor * k0
}

void Linear_Convection_1d::lax_wendroff()
{
    solution_new[0] = solution_old[0] - 0.5 * cfl * (solution_old[1] - solution_old[n_points - 1]) + 0.5 * cfl * cfl * (solution_old[n_points - 1] - 2 * solution_old[0] + solution_old[1]);
    for (int j = 1; j < n_points - 1; ++j) //Loop over grid points
    {
        solution_new[j] = solution_old[j] - 0.5 * cfl * (solution_old[j + 1] - solution_old[j - 1]) + 0.5 * cfl * cfl * (solution_old[j - 1] - 2 * solution_old[j] + solution_old[j + 1]);
    }
    solution_new[n_points - 1] = solution_old[n_points - 1] - 0.5 * cfl * (solution_old[0] - solution_old[n_points - 2]) + 0.5 * cfl * cfl * (solution_old[n_points - 2] - 2 * solution_old[n_points - 1] + solution_old[0]);
};

//This computes solution_new = u^{n+1} = u^n + dt/6 * (k1 + 2.0*k2 + 2.0*k3 + k4)
void Linear_Convection_1d::compute_solution_new_using_ki()
{
    for (int i = 0; i < n_points; i++)
        solution_new[i] = solution_old[i] + dt / 6 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
}

void Linear_Convection_1d::rk4_solver()
{
    //This innocent step was what was causing trouble
    /* solution_old = solution_new; */

    //k1 = rhs_function(solution_old)
    rhs_function(solution_old, k1);
    //Temporarily putting u^{n+1} = solution_new = solution_old + dt/2 * k1
    temporary_update_solution(dt / 2, k1, solution_new);
    //So, computing k2 = rhs_function(u^n + dt/2 * k1)
    rhs_function(solution_new, k2);
    //Similarly, temporarily putting u^{n+1}=solution_new = u^n + dt/2 *k2
    temporary_update_solution(dt / 2, k2, solution_new);
    //Computing k3 = rhs_function(solution_old + dt/2 *k2)
    rhs_function(solution_new, k3);
    //Temporarily putting u^{n+1} = solution_new = u^n + dt * k3
    temporary_update_solution(dt, k3, solution_new);
    //Computing k4 = rhs_function(solution_old + dt *k3)
    rhs_function(solution_new, k4);

    //This computes solution_new = u^{n+1} = u^n + dt/6 * (k1 + 2.0*k2 + 2.0*k3 + k4)
    compute_solution_new_using_ki();
};

void Linear_Convection_1d::output_error_and_results(int time_step_number)
{
    if (time_step_number == 0)
        solution_exact = initial_data;
    else
    {
        switch (initial_data_indicator)
        {
        case 0:
            for (int j = 0; j < n_points; j++)
            {
                solution_exact[j] = sin(2 * M_PI * (grid[j] - coefficient * dt * time_step_number));
            }
            break;
        case 1:
            for (int j = 0; j < n_points; j++)
            {
                solution_exact[j] = hat_function(grid[j] - coefficient * dt * time_step_number);
            }
            break;
        case 2:
            for (int j = 0; j < n_points; j++)
            {
                solution_exact[j] = step_function(grid[j] - coefficient * dt * time_step_number);
            }
            break;
        default:
            cout << "You entered the wrong initial_data_indicator ";
            assert(false);
        }
    }
    string solution_file_name = "solution_", error_file_name = "error_";
    solution_file_name += to_string(time_step_number), error_file_name += to_string(time_step_number);
    solution_file_name += "_" + method, error_file_name += "_" + method;
    solution_file_name += ".txt", error_file_name += ".txt";
    ofstream output_solution, output_error;
    output_solution.open(solution_file_name), output_error.open(error_file_name);
    for (int j = 0; j < n_points; j++)
    {
        error[j] = abs(solution_old[j] - solution_exact[j]);
        if (j > 0)
            output_solution << "\n", output_error << "\n";
        output_solution << grid[j] << " " << solution_old[j], output_error << grid[j] << " " << error[j];
    }
    output_solution.close();
}

void Linear_Convection_1d::run_and_output_results()
{
    make_grid();
    set_initial_data();
    solution_old = initial_data;
    int time_step_number = 0; //The bug I had made here was that I was running the solver
    //even for time_step_number = 0, which is wrong.
    output_error_and_results(time_step_number);
    while (time_step_number * dt < running_time)
    {
        time_step_number += 1;
        if (method == "rk4")
            rk4_solver();
        else if (method == "lw")
            lax_wendroff();
        else
            assert(false);
        solution_old = solution_new;
        output_error_and_results(time_step_number);
    }
    cout << "We made " << time_step_number << " steps.";
}

int main()
{
    string method;
    double n_points = 60;
    double cfl;
    double running_time;
    int initial_data_indicator;
    cout << "Please type 'lw' for Lax-Wendroff and 'rk4' for Runge-Kutta 4." << endl;
    cin >> method;
    cout << "Please enter the cfl number " << endl;
    cin >> cfl;
    cout << "How much t would you like the solution to run?" << endl;
    cin >> running_time;
    cout << "Please pick where you'd like sin(x), hat or step function as initial data by typing 0,1 or 2 respectively " << endl;
    cin >> initial_data_indicator;
    Linear_Convection_1d solver(n_points, cfl, method, running_time, initial_data_indicator);
    solver.run_and_output_results();
}


double Linear_Convection_1d::step_function(double grid_point)
{
    double value;
    if (grid_point < x_min + (x_max - x_min) / 4.0 || grid_point > x_max - (x_max - x_min) / 4.0)
        value = 0;
    else
        value = 1.0;
    return value;
}