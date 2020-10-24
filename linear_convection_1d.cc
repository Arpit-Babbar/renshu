#include <iostream>
#include <fstream>
#include <cmath>
#include <vector> 

using namespace std;

class Linear_Convection
{
    public:
      Linear_Convection(double x_min0, double x_max0, int n_points0, double cfl0, double coefficient0,vector<double> initial_data);
      void run();
    
    private:
      void make_grid();
      //void setup_system(); Not needed
      void solve_system(); //This should give the solution at time t = n_iterations * time_step
      
      double coefficient;double x_min;double x_max; //PDE parameters

      double n_points; double cfl; double time_step; int n_iterations; //Numerical algorithm parameters
      vector<double> grid;

      vector<double> solution_old; //solution at time t-1
      vector<double> solution_t; //solution at time t
};
Linear_Convection::Linear_Convection(double x_min0, double x_max0, int n_points0, double cfl0,double coefficient0, vector<double> initial_data)
  : x_min(x_min0), x_max(x_max0), n_points(n_points0) , cfl(cfl0)
   , coefficient(coefficient0)
   , time_step(n_points/(coefficient * n_points))
   , solution_t(n_points) //specifying size
   , solution_old(initial_data) //specifying the full vector
{}
void Linear_Convection::make_grid()
{
    for (int j = 0; j < n_points ; j++)
    {
        double h = (x_max - x_min)/n_points;
        grid.push_back(x_min + j*h);
    }   
    for (int j = 0; j < grid.size() ; j++)
    {
        cout << grid.at(j) << endl;
    }

}



void Linear_Convection::solve_system()
{
    for (int j = 0; j < n_points; ++j) //Loop over grid points
    {
        solution_t[j] = solution_old[j] - 0.5 * cfl * (solution_old[j+1] - solution_old[j-1]) + 0.5 * cfl * cfl * (solution_old[j-1] - 2 * solution_old[j] + solution_old[j+1]);
        
    }
    solution_old = solution_t;
}

void Linear_Convection::run()
{
    make_grid();
    n_iterations = 10;
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
    //set parameters
    double x_min = 0; double x_max = 1.0; double coefficient = 1.0;
    double n_points = 30; double cfl = 0.5;
    vector<double> initial_data(n_points);
    for (auto value_i = initial_data.begin(); value_i != initial_data.end(); ++value_i) 
    {
        *value_i = 1.0;     
    }
    //make_grid
    Linear_Convection solver(x_min,x_max,n_points,cfl, coefficient, initial_data);
    solver.run();
    return 0;
}
