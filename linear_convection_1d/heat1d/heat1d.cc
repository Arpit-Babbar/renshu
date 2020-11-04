//u_t + u_xx = 0, u(x,0) = sin(2*pi*x) with periodic and zero bc.
//exact solution is exp(-a 4 pi^2 t)sin(2*pi*x)

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <string>

#include <algorithm>
#include <functional> //Used to define addition of vectors

using namespace std;

double max_element(vector<double> &v)
{
    double a = 0.0;
    for (int i = 0; i < v.size(); i++)
        a = max(a, v[i]);
    return a;
}

class Linear_Convection_1d
{
public:
    Linear_Convection_1d(const double n_points, const double cfl, string method, const double running_time, int initial_data_indicator); //input parameters
    //and which method to use - Lax-Wendroff or RK4
    void run(); //true when output is to be given, and false when it doesn't;
    void evaluate_error_and_output(const int time_step_number);
    double get_max_error()
    {
        return max_element(error);
    }
    double get_l2_error()
    {
        double a = 0.0;
        for (int j = 0; j < n_points; j++)
        {
            a = a + error[j] * error[j] * h;
        }
        return sqrt(a);
    }
    void indicate_to_output() { output_indicator = 1; }

private:
    void make_grid(); // Grid is needed for writing output, and defining exact solution.
    void set_initial_data();
    void temporary_update_solution(const double factor, const vector<double> &u, vector<double> &k);
    //We'd use this to compute the rk4 slopes k_i's and temporarily store them in solution_new.
    //More precisely, it does k = solution_old + factor * u
    void temporary_update_solution(const double factor0, const vector<double> &u0,
                                   const double factor1, const vector<double> &u1, vector<double> &k);
    void rk4_solver(); //Gives the solution at next time step using RK4
    void rk3_solver();
    void rk2_solver();
    void ftcs();                                                   //Gives the solutions at next time step using Lax-Wendroff
    void rhs_function(const vector<double> &u, vector<double> &k); //This gives RHS of the system of ODEs

    double coefficient = 1.0;

    double x_min = 0.0;
    double x_max = 1.0;

    double interval_part(double); //If a point is not in the interval [x_min,x_max], this function sends it there
    //by translating in x_max-x_min sized steps.

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
    void compute_solution_new_using_ki(int);

    string method;
    int initial_data_indicator; //Kept in case we want to add multiple initial datas
    int output_indicator = 0;
};

Linear_Convection_1d::Linear_Convection_1d(double n_points, double cfl, string method, double running_time, int initial_data_indicator) : n_points(n_points), cfl(cfl), running_time(running_time), method(method), initial_data_indicator(initial_data_indicator)
{
    h = (x_max - x_min) / (n_points - 1);
    dt = cfl * h * h / coefficient;
    std::cout << "Spatial grid points are at gap h = " << h << endl;
    std::cout << "Time step dt = " << dt << endl;
    grid.resize(n_points);
    initial_data.resize(n_points);
    error.resize(n_points);
    solution_old.resize(n_points);
    solution_new.resize(n_points);
    solution_exact.resize(n_points);
    std::cout << "Number of spatial grid points is " << n_points << endl;
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

void Linear_Convection_1d::set_initial_data()
{
    for (int i = 0; i < n_points; i++)
        initial_data[i] = sin(2 * M_PI * grid[i] / (x_max - x_min));
}

//This computes the rhs of the system of ODEs on which we apply rk4.
void Linear_Convection_1d::rhs_function(const vector<double> &u, vector<double> &k)
{
    if (u.size() != k.size() || k.size() != solution_old.size())
    {
        std::cout << "You have used rhs_function to do rhs_funciton(&u,&k)"
             << " with inappropriately sized u,k. Size of your u is " << u.size() << ", size of k is " << k.size()
             << ". Both the sizes should equal " << solution_old.size() << endl;
        assert(false);
    }
    k[0] = (u[n_points - 1] - 2 * u[0] + u[1]) / (h * h); //left end point
    for (int j = 1; j < n_points - 1; j++)
    {
        k[j] = (u[j - 1] - 2 * u[j] + u[j + 1]) / (h * h);
    }
    k[n_points - 1] = (u[n_points - 2] - 2 * u[n_points - 1] + u[0]) / (h * h); //right end point
};

//k = solution_old + factor*u
void Linear_Convection_1d::temporary_update_solution(const double factor, const vector<double> &k0, vector<double> &k)
{
    if (k0.size() != k.size() || k.size() != solution_old.size())
    {
        std::cout << "You have used temporary_update_solution to do 'k = solution_old + factor * k0 "
             << "with inappropriately sized k0,k. Size of your k0 is " << k0.size() << ", size of k is " << k.size()
             << ". Both the sizes should equal " << solution_old.size() << endl;
        assert(false);
    }
    for (int i = 0; i < k.size(); i++)
    {
        k[i] = solution_old[i] + factor * k0[i];
    } //k = solution_old + factor * k0
}
void Linear_Convection_1d::temporary_update_solution(const double factor0, const vector<double> &k0, const double factor1, const vector<double> &k1, vector<double> &k)
{
    if (k0.size() != k.size() || k.size() != solution_old.size())
    {
        std::cout << "You have used temporary_update_solution to do 'k = solution_old + factor * k0 "
             << "with inappropriately sized k0,k. Size of your k0 is " << k0.size() << ", size of k is " << k.size()
             << ". Both the sizes should equal " << solution_old.size() << endl;
        assert(false);
    }
    for (int i = 0; i < k.size(); i++)
    {
        k[i] = solution_old[i] + factor0 * k0[i] + factor1 * k1[i];
    } //k = solution_old + factor * k0
}

void Linear_Convection_1d::ftcs()
{
    solution_new[0] = 0.0;
    for (int j = 1; j < n_points - 1; ++j) //Loop over grid points
    {
        solution_new[j] = cfl * solution_old[j - 1] + (1 - 2 * cfl) * solution_old[j] + cfl * solution_old[j + 1];
    }
    solution_new[n_points - 1] = 0.0;
};

//This computes solution_new = u^{n+1} = u^n + dt/6 * (k1 + 2.0*k2 + 2.0*k3 + k4)
void Linear_Convection_1d::compute_solution_new_using_ki(int order)
{
    if (order == 4)
    {
        for (int i = 0; i < n_points; i++)
            solution_new[i] = solution_old[i] + dt / 6 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }
    if (order == 3)
    {
        for (int i = 0; i < n_points; i++)
            solution_new[i] = solution_old[i] + dt / 6 * (k1[i] + 4.0 * k2[i] + k3[i]);
    }
    if (order == 2)
    {
        for (int i = 0; i < n_points; i++)
            solution_new[i] = solution_old[i] + dt / 2 * (k1[i] + k2[i]);
    }
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
    compute_solution_new_using_ki(4.0);
};

void Linear_Convection_1d::rk3_solver()
{
    //This innocent step was what was causing trouble
    /* solution_old = solution_new; */

    //k1 = rhs_function(solution_old)
    rhs_function(solution_old, k1);
    //Temporarily putting u^{n+1} = solution_new = solution_old + dt/2 * k1
    temporary_update_solution(dt / 2, k1, solution_new);
    //So, computing k2 = rhs_function(u^n + dt/2 * k1)
    rhs_function(solution_new, k2);
    //Similarly, temporarily putting u^{n+1}=solution_new = u^n -k1 * dt + 2 k2 * dt
    temporary_update_solution(-dt, k1, 2 * dt, k2, solution_new);
    //Computing k3 = rhs_function(solution_old - dt k1 + 2 dt k2)
    rhs_function(solution_new, k3);

    //This computes solution_new = u^{n+1} = u^n + dt/6 * (k1 + 2.0*k2 + 2.0*k3 + k4)
    compute_solution_new_using_ki(3.0);
};

void Linear_Convection_1d::rk2_solver()
{
    //This innocent step was what was causing trouble
    /* solution_old = solution_new; */

    //k1 = rhs_function(solution_old)
    rhs_function(solution_old, k1);
    //Temporarily putting u^{n+1} = solution_new = solution_old + dt * k1
    temporary_update_solution(dt, k1, solution_new);
    //So, computing k2 = rhs_function(u^n + dt * k1)
    rhs_function(solution_new, k2);

    //This computes solution_new = u^{n+1} = u^n + dt/6 * (k1 + 2.0*k2 + 2.0*k3 + k4)
    compute_solution_new_using_ki(2.0);
};

void Linear_Convection_1d::evaluate_error_and_output(int time_step_number)
{
    if (time_step_number == 0)
        solution_exact = initial_data;
    else
    {
        for (int j = 0; j < n_points; j++)
        {
            solution_exact[j] = exp(-4 * coefficient * M_PI * M_PI * dt * time_step_number) * sin(2 * M_PI * (grid[j] - dt * time_step_number) / (x_max - x_min));
        }

        for (int j = 0; j < n_points; j++)
        {
            error[j] = abs(solution_old[j] - solution_exact[j]);
        }
        if (output_indicator == true) //Give output only if asked for
        {
            string solution_file_name = "solution_", error_file_name = "error_";
            solution_file_name += to_string(time_step_number), error_file_name += to_string(time_step_number);
            solution_file_name += "_" + method, error_file_name += "_" + method;
            solution_file_name += ".txt", error_file_name += ".txt";
            ofstream output_solution, output_error;
            output_solution.open(solution_file_name), output_error.open(error_file_name);
            for (int j = 0; j < n_points; j++)
            {
                output_solution << grid[j] << " " << solution_old[j] << " " << solution_exact[j] << "\n";
                output_error << grid[j] << " " << error[j] << "\n";
            }
            output_solution.close();
        }
    }
}

    void Linear_Convection_1d::run()
    {
        make_grid();
        set_initial_data();
        solution_old = initial_data;
        int time_step_number = 0; //The bug I had made here was that I was running the solver
        //even for time_step_number = 0, which is wrong.
        evaluate_error_and_output(time_step_number);
        while (time_step_number * dt < running_time)
        {
            time_step_number += 1;
            if (method == "rk4")
                rk4_solver();
            else if (method == "rk3")
                rk3_solver();
            else if (method == "rk2")
                rk2_solver();
            else if (method == "ftcs")
                ftcs();
            else
                assert(false);
            solution_old = solution_new;
            evaluate_error_and_output(time_step_number);
        }
        cout << "We made " << time_step_number << " steps.";
    }

    void run_and_get_output(double n_points, double cfl, string method, double running_time, int initial_data_indicator, double tolerance)
    {
        ofstream error_rate;
        error_rate.open("error_rate.txt");
        Linear_Convection_1d solver(n_points, cfl, method, running_time, initial_data_indicator);
        vector<double> errors;
        double iteration_number = 0;
        solver.run();
        cout << solver.get_max_error();
        while (solver.get_max_error() > tolerance)
        {
            errors[iteration_number] = solver.get_max_error();
            error_rate << 1 / n_points << " " << errors[iteration_number] << "\n";
            n_points = 2.0 * n_points;
            solver = Linear_Convection_1d(n_points, cfl, method, running_time, initial_data_indicator);
            if (iteration_number > 0)
            {
                std::cout << "Rate of convergence checked at iteration number " << iteration_number + 1;
                std::cout << " is " << abs(log(errors[iteration_number] / errors[iteration_number - 1])) / log(2.0) << endl;
            }
            solver.run(); //running without generating output
            iteration_number++;
        }
        error_rate.close();
    solver.indicate_to_output();
    solver.run();
    cout << "It took "<< iteration_number << " iterations to get the Linfty error below " << tolerance<< " and precisely at "<< solver.get_max_error()<< endl << endl;
    }
    

    int main(int argc, char **argv)
    {
        if (argc != 5)
        {
            std::cout << "Incorrect arguments. Kindly give the arguments in the following format. " << endl;
            std::cout << "./output ftcs/rk4/rk3/rk2(for the respective method) cfl running_time error_tolerance" << endl;
            assert(false);
        }
        string method = argv[1];
        std::cout << "You have entered the method to be " << method << endl;
        double n_points = 5.0;
        double cfl = stod(argv[2]);
        std::cout << "You have entered the cfl number to be " << cfl << endl;
        double running_time = stod(argv[3]);
        std::cout << "You have picked the running_time to be " << running_time << endl;
        int initial_data_indicator = 0;
        double tolerance = stod(argv[4]);
        cout << "You have entered the tolerance to be "<< tolerance;
        run_and_get_output(n_points, cfl, method, running_time, initial_data_indicator, tolerance);
    }
