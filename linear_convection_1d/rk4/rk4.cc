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
    std::vector<T> result;
    result.reserve(V.size());

    for (int i = 0; i < V.size(); i ++)
    {
        result[i] = a*V[i];
    }
    return result;
}

vector<double> scalar_multiply(double scalar,vector<double> u)
{
    vector<double> result(u.size());
    for (int i = 0; i < u.size(); i++)
    {
        result[i] = scalar * u[i];
    };
    return result;
}
class advection_rk4
{
    public:
    advection_rk4();
    void run();

    private:
    void solve_system(); //Gives the solution at next time step;
    void make_grid();
    vector<double> f(double, vector<double>); //This gives RHS of the system of ODEs
    //on which we apply Runge-Kutta.
    double coefficient = 1.0;

    double x_min = 0.0; double x_max = 1.0;

    vector<double> grid;

    vector<double> solution_old;
    vector<double> solution_new;

    double n_points; double h; double dt; double cfl;
};

advection_rk4::advection_rk4():n_points(n_points)
            , grid(n_points), solution_old(n_points), solution_new(n_points)
            , h((x_max - x_min)/n_points),dt(dt)
            , cfl(abs(coefficient)* dt /h)
{};

void advection_rk4::make_grid()
{
    for (int i = 0; i < n_points; i++)
    {
        grid[i] = x_min + i * h;
    };
};

vector<double> advection_rk4::f(double h, vector<double> u)
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

void advection_rk4::solve_system()
{
    vector<double> initial_data(n_points);

    for (int i = 0; i< n_points; i++) //Setting initial data to be constant 1
    {
        initial_data[i] = 1.0;
    }

    for (int j =0 ; j<n_points; j++)
    {
        if (j==0)
        {
            solution_new = initial_data;
            solution_old = solution_new;
            break;
        }
        vector<double> k1; vector<double> k2; vector<double> k3;vector<double> k4;
        k1 = f(h,solution_old);k2 = f(h, solution_old + dt/2 * k1);
        k3 = f(h,solution_old + dt/2 * k2);
        k4 = f(h,solution_old + dt * k3);

        solution_new = solution_old;
        vector<double> values_to_be_added(n_points); //This will contain all values we add
        //to progress solution by one time step. This step may need to be removed
        //for optimization.
        solution_new = solution_new + dt/6 * (k1 + 2.0*k2 + 2.0*k3 + k4);
    };
};

void advection_rk4::run()
{
    make_grid();
    int n_iterations = 10;
    for (int i = 0; i < n_iterations; i++)
    {
        if (i == 0)
        {
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
            solve_system();
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
    double n_points = 30; double dt = 0.01; //This makes cfl = 0.3
    double a = 3; vector<double> b = {2};
    b = a*b;
    std::cout <<b[0];
}