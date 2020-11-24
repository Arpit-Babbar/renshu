#include <vector>
#include <iostream>
#include <fstream>

#include <cassert>
using namespace std;
///Defining y = x1 + a2*x2
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

//Defining y = a1*x1 + a2*x2
void add(const double a1, const vector<double> &x1,
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
        y[i] = a1*x1[i] + a2*x2[i];
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
