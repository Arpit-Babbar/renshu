#include <iostream>
#include <fstream>
#include <cmath>
#include <vector> 
#include <cassert>

#include<algorithm>
#include<functional> //Used to define addition of vectors

using namespace std;

//Got vector adding from here https://stackoverflow.com/questions/3376124/how-to-add-element-by-element-of-two-stl-vectors
//Used it to make scalar multiplication
template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
    if (!(a.size() == b.size()))
    {
        cout << "These are the unequal vectors "<< endl;
        for (int i = 0; i < a.size(); i++)
            {cout <<  a[i] << " ";}
        cout << endl;
        cout << "Size of b is " << b.size() << ", and entries of b are "<< endl;
        for (int i = 0; i < b.size(); i++)
            {cout << b[i] << " ";}
        cout << endl;
        
    }
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::plus<T>());
    return result;
}

template <typename T>
std::vector<T> operator*(const T a, const std::vector<T>& V) //scalar multiplication
{
    std::vector<T> result(V.size());

    for (int i = 0; i < V.size(); i ++)
    {
        result[i] = a*V[i];
    }
    return result;
}

class Linear_Convection_1d
{
    public:
    Linear_Convection_1d(double n_points, double dt,string method);
    void run();

    private:
    void make_grid();
    void rk4_solver(); //Gives the solution at next time step;
    void lax_wendroff();
    vector<double> f(double, vector<double>); //This gives RHS of the system of ODEs
    //on which we apply Runge-Kutta.
    double coefficient = 1.0;

    double x_min = 0.0; double x_max = 1.0;

    vector<double> grid;

    vector<double> solution_old;
    vector<double> solution_new;

    double n_points; double h; double dt; double cfl;

    string method;
};

Linear_Convection_1d::Linear_Convection_1d(double n_points, double dt, string method):n_points(n_points)
            , grid(n_points), solution_old(n_points), solution_new(n_points)
            , h((x_max - x_min)/n_points),dt(dt)
            , cfl(abs(coefficient)* dt /h)
            , method(method)
{};

void Linear_Convection_1d::make_grid()
{
    for (int i = 0; i < n_points; i++)
    {
        grid[i] = x_min + i * h;
    };
};

vector<double> Linear_Convection_1d::f(double h, vector<double> u)
{
    vector<double> right_hand_side(u.size());
    for (int j = 0; j < n_points; j++)
    {
        if (j==0)
        right_hand_side[j] = -coefficient * (u[j+1] - u[n_points - 1]) / (2 * h);
        else if (j > 0 && j < n_points - 1)
        right_hand_side[j] = -coefficient * (u[j+1] - u[j-1]) / (2 * h);
        else if (j == n_points - 1)
        right_hand_side[j] = -coefficient * (u[0] - u[j-1]) / (2 * h);
    }
    return right_hand_side;
};
void Linear_Convection_1d::lax_wendroff()
{
    for (int j = 0; j < n_points; ++j) //Loop over grid points
    {
        if (j== 0)  //Periodic boundary condition
        solution_new[j] = solution_old[j]
                        - 0.5 * cfl * (solution_old[j+1] - solution_old[n_points - 1])
                        + 0.5 * cfl * cfl * (solution_old[n_points-1] - 2 * solution_old[j] + solution_old[j+1]);
        else if (j > 0 && j < n_points - 1)
        solution_new[j] = solution_old[j]
                        - 0.5 * cfl * (solution_old[j+1] - solution_old[j-1])
                        + 0.5 * cfl * cfl * (solution_old[j-1] - 2 * solution_old[j] + solution_old[j+1]);
        else if (j == n_points - 1) //periodic boundary condition
        solution_new[j] = solution_old[j]
                        - 0.5 * cfl * (solution_old[0] - solution_old[j-1])
                        + 0.5 * cfl * cfl * (solution_old[j-1] - 2 * solution_old[j] + solution_old[0]);
                
    }
    solution_old = solution_new;
};
void Linear_Convection_1d::rk4_solver()
{
        vector<double> k1(n_points); vector<double> k2(n_points); vector<double> k3(n_points);vector<double> k4(n_points);
        k1 = f(h,solution_old);
        vector<double> alpha = dt/2 * k1;
        k2 = f(h, solution_old + dt/2 * k1);
        k3 = f(h,solution_old + dt/2 * k2);
        k4 = f(h,solution_old + dt * k3);

        solution_new = solution_old;
        solution_new = solution_new + dt/6 * (k1 + 2.0*k2 + 2.0*k3 + k4);
        solution_old = solution_new;
};

void Linear_Convection_1d::run()
{
    make_grid();
    vector<double> initial_data(n_points);
    for (int i = 0; i< n_points; i++) //Setting initial data to be constant 1
    {
        initial_data[i] = i*h;
    }
    int n_iterations = 10;
    for (int i = 0; i < n_iterations; i++)
    {
        if (i == 0)
        {
            solution_old = initial_data;
            string file_name = "initial_solution.txt";
            ofstream output_solution;
            output_solution.open (file_name);
            for (int j = 0; j < n_points; j++)
            {
            if (j>0)
            output_solution << "\n";
            output_solution << grid[j] << " " << solution_old[j];
            }
            output_solution.close();
        }
        else
        {
            if (method == "rk4")
            rk4_solver();
            else if (method == "lw")
            lax_wendroff();
            else
            assert(false);
            string file_name = "solution_";
            file_name += to_string(i); file_name += ".txt";
            ofstream output_solution;
            output_solution.open (file_name);
            for (int j = 0; j < n_points; j++)
            {
            if (j>0)
            output_solution << "\n";
            output_solution << grid[j] << " " << solution_old[j];
            }
            output_solution.close();
        }
        
    }
}

int main()
{
    string method;
    double n_points = 30; double dt = 0.01; //This makes cfl = 0.3
    cout << "Please type 'lw' for Lax-Wendroff and 'rk4' for Runge-Kutta 4."<<endl;
    cin >> method;
    Linear_Convection_1d solver(n_points, dt,method);
    solver.run();
}