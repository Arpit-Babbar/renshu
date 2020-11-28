#include <vector>
#include <iostream>
#include <fstream>

#include <cassert>
using namespace std;
///Defining y = x1 + a2*x2
void add(const vector<double> &x1,
         const double a2, const vector<double> &x2,
         vector<double> &y);

//Defining y = a1*x1 + a2*x2
void add(const double a1, const vector<double> &x1,
         const double a2, const vector<double> &x2,
         vector<double> &y);

//Defining y = x1 + a2*x2 + a3*x3
void add(const vector<double> &x1,
         const double a2, const vector<double> &x2,
         const double a3, const vector<double> &x3,
         vector<double> &y);

//output vectors as file_name.txt with columns in following format
//vector1 vector2 vector3
void output_vectors_to_file(string file_name, vector<double> &grid,
                            vector<double> &solution_old, vector<double> &solution_exact);

//output vectors as file_name.txt with columns in following format
//vector1 vector2
void output_vectors_to_file(string file_name, vector<double> &grid, vector<double> &solution_old);