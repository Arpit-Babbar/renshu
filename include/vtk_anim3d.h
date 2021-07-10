#ifndef __VTK_ANIM_H__
#define __VTK_ANIM_H__
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cassert>

#include "array3d.h"

using namespace std;

// Ex: get_filename("sol_",3,42) should return "sol_042.vtk"
string get_filename(const string base_name,
                    const int ndigits,
                    const int c);
//Forms grid by taking tensor product of grid_x and grid_y
//Solution is defined on the respective grid points.
void write_rectilinear_grid(vector<double> &grid_x,
                            vector<double> &grid_y,
                            vector<double> &grid_z,
                            Array3D &solution,
                            double t,
                            int c, //Cycle number
                            string filename);

void vtk_anim_sol(vector<double> &grid_x,vector<double> &grid_y,
                  Array3D& solution,
                   double t,
                  int time_step_number,
                  string filename);

//Double input

void write_rectilinear_grid(vector<double> &grid_x,
                            vector<double> &grid_y,
                            vector<double> &grid_z,
                            Array3D& solution, Array3D& solution_exact,
                            double t,
                            int c, //Cycle number
                            string filename);

void vtk_anim_sol(vector<double> &grid_x,vector<double> &grid_y,
                  Array3D& solution, Array3D& solution_exact,
                   double t,
                  int time_step_number,
                  string filename);
#endif
