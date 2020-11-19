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
    void rk4_solver(); //Gives the solution at next time step using RK4
    void rk3_solver();
    void rk3_solver_new();
    void rk2_solver();
    void lax_wendroff();                                           
    void rhs_function(); //This gives RHS of the system of ODEs and stores it to where rhs points.
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

    vector<double> *rhs;

    vector<double> solution_exact; //Exact solution at present time step

    vector<double> error; //This will store the maximum error at a grid point in all time-steps
    double l2_error;

    //h denotes the spatial distance between grid points
    double n_points, h, dt,t, cfl, running_time; //h = 1/n_points just included for easy typing

    //slopes needed by rk4
    vector<double> temp;

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
        default:
            cout << "You entered the wrong initial_data_indicator ";
            assert(false);
        }
    }
}

//This computes the rhs of the system of ODEs on which we apply rk4.
void Linear_Convection_1d::rhs_function()
{
    (*rhs)[0] = -(solution[1] - solution[n_points - 1]) / (2.0 * h); //left end point
    for (int j = 1; j < n_points - 1; j++)
    {
        (*rhs)[j] = -(solution[j + 1] - solution[j - 1]) / (2.0 * h);
    }
    (*rhs)[n_points - 1] = -(solution[0] - solution[n_points - 2]) / (2.0 * h); //right end point
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
    //Since rhs_function gives its output to (*rhs), and we are interested
    //in updating temp, we shall do this sneaky thing.
    rhs = &temp;
    //We know the rk4 time-stepping formula explicitly. We only need temp
    //because steps like solution = rhs_function(solution) do not work 
    //because of the way we have defined rhs_function.
    rhs_function();//temp = rhs(solution) = rhs(solution_old)
    
    add(solution_old,dt/4.0,temp,solution); 
    //Temporarily putting u^{n+1} = u^n + dt/4*temp

    rhs_function();//temp = rhs(solution)
    add(solution_old,dt/3.0,temp,solution);
    //u^{n+1} = u^n + dt/3 * temp

    rhs_function(); //temp = rhs(solution)
    add(solution_old,dt/2.0,temp,solution);

    rhs_function();
    add(solution_old,dt,temp,solution);
}

void Linear_Convection_1d::rk3_solver()
{
    rhs = &temp;
    //k1 = rhs_function(solution) = rhs_function(solution_old)
    rhs_function();
    //Temporarily putting u^{n+1} = solution_new = solution_old + dt/2 * k1
    add(solution_old, dt / 3.0, temp, solution);
    //So, computing k2 = rhs_function(u^n + dt/2 * k1)
    rhs_function();
    add(solution_old, 0.5 * dt, temp, solution);
    rhs_function();
    add(solution_old, dt, temp, solution);
}

void Linear_Convection_1d::rk3_solver_new()
{
    rhs = &temp;
    //k1 = rhs_function(solution) = rhs_function(solution_old)
    rhs_function();
    //Temporarily putting u^{n+1} = solution_new = solution_old + dt/2 * k1
    add(solution_old, dt / 6.0, temp, solution);
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
    rhs_function();
    //Temporarily putting u^{n+1} = solution = solution_old + 0.5 * dt * k1
    add(solution_old, 0.5 * dt, temp, solution);
    //Next, putting k1 = rhs_function(solution) = rhs_function(solution+old + 0.5)
    rhs_function();
    add(solution_old,dt, temp, solution);
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
    set_initial_solution(); //sets solution to be the initial data
    int time_step_number = 0; //The bug I had made here was that I was running the solver
    //even for time_step_number = 0, which is wrong.
    evaluate_error_and_output_solution(time_step_number);
    while (t < running_time) //compute solution at next time step using solution_old
    {
        solution_old = solution;//update solution_old to be used at next time step
        if (method == "rk4")
            rk4_solver();
        else if (method == "rk3")
            rk3_solver();
        else if (method == "rk3_new")
            rk3_solver_new();
        else if (method == "rk2")
            rk2_solver();
        else if (method == "lw")
            lax_wendroff();
        else
            assert(false);
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
        if (refinement_level == max_refinements + 1) solver.output_final_error(); //I couldn't figure a
        //way to keep it outside the while loop.
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
        cout << "./output lw/rk4/rk3/rk3_new/rk2(for the respective method) cfl running_time 0/1/2(for " ;
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
    unsigned int max_refinements = stoi(argv[5]);
    cout << "max_refinements = " << max_refinements <<endl;
    run_and_get_output(n_points, cfl, method, running_time, initial_data_indicator, max_refinements);
}