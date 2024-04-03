#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cassert>

#include "array3d.h"
#include "vtk_anim3d.h"

using namespace std;

// Ex: get_filename("sol_",3,42) should return "sol_042.vtk"
string get_filename(const string base_name,
                    const int ndigits,
                    const int c)
{
   if(c > pow(10,ndigits)-1)
   {
      cout << "get_filename: Not enough digits !!!\n";
      cout << "ndigits= " << ndigits << endl;
      cout << "c      = " << c << endl;
      exit(0);
   }

   string name = base_name;
   // first pad with zeros
   int d = 1;
   if(c > 0) d = int(floor(log10(c))) + 1;
   for(int i=0; i<ndigits-d; ++i)
      name += "0";
   name += to_string(c) + ".vtk";
   return name;
}
//Forms grid by taking tensor product of grid_x and grid_y
//Solution is defined on the respective grid points.
void write_rectilinear_grid(vector<double> &grid_x,
                            vector<double> &grid_y,
                            vector<double> &grid_z,
                            Array3D &solution,
                            double t,
                            int c, //Cycle number
                            string filename)
{
   const int nx = solution.sizex();
   const int ny = solution.sizey();
   const int nz = solution.sizez(); // We have a 2d grid
   ofstream fout;
   fout.open(filename);
   fout << "# vtk DataFile Version 3.0" << endl;
   fout << "Cartesian grid" << endl;
   fout << "ASCII" << endl;
   fout << "DATASET RECTILINEAR_GRID" << endl;
   fout << "FIELD FieldData 2" << endl;
   fout << "TIME 1 1 double" << endl;
   fout << t << endl;
   fout << "CYCLE 1 1 int" << endl;
   fout << c << endl;
   fout << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
   fout << "X_COORDINATES " << nx << " float" << endl;
   for(int i=0; i<nx; ++i)
      fout << grid_x[i] << " ";
   fout << endl;
   fout << "Y_COORDINATES " << ny << " float" << endl;
   for(int j=0; j<ny; ++j)
      fout << grid_y[j] << " ";
   fout << endl;
   fout << "Z_COORDINATES " << nz << " float" << endl;
   for(int k=0; k<nz; ++k)
      fout << grid_z[k] << " ";
   fout << endl;

   fout << "POINT_DATA " << nx*ny*nz << endl;
   fout << "SCALARS density float" << endl;
   fout << "LOOKUP_TABLE default" << endl;
   // no need for k-loop since nk=1
   for(int i=0; i<nx; ++i)
   {
      for(int j=0; j<ny; ++j)
      {
         for(int k=0; k<nz; ++k)
            fout << solution(i,j,k) << " ";
         fout << endl;
      }
   }
}
void vtk_anim_sol(vector<double> &grid_x,vector<double> &grid_y,
                  vector<double> &grid_z,
                  Array3D& solution,
                   double t,
                  int time_step_number,
                  string filename)
{
  filename = filename+"_";
  filename = get_filename(filename,3,time_step_number);
  write_rectilinear_grid(grid_x, grid_y, grid_z, solution, t, time_step_number,
                         filename);
}


// Double input

void write_rectilinear_grid(vector<double> &grid_x,
                            vector<double> &grid_y,
                            vector<double> &grid_z,
                            Array3D& solution, Array3D& solution_exact,
                            double t,
                            int c, //Cycle number
                            string filename)
{
   const int nx = solution.sizex();
   const int ny = solution.sizey();
   const int nz = solution.sizez(); // We have a 2d grid
  /*  fout.open(filename)
      write_grid(fout,nx,ny,dx,dy)
      write_sol(fout,nxny,solution,"solution")
      write_sol(fout,nx,ny,solution_exact,"solution_exact")
  */
   ofstream fout;
   fout.open(filename);
   fout << "# vtk DataFile Version 3.0" << endl;
   fout << "Cartesian grid" << endl;
   fout << "ASCII" << endl;
   fout << "DATASET RECTILINEAR_GRID" << endl;
   fout << "FIELD FieldData 2" << endl;
   fout << "TIME 1 1 double" << endl;
   fout << t << endl;
   fout << "CYCLE 1 1 int" << endl;
   fout << c << endl;
   fout << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
   fout << "X_COORDINATES " << nx << " float" << endl;
   for(int i=0; i<nx; ++i)
      fout << grid_x[i] << " ";
   fout << endl;
   fout << "Y_COORDINATES " << ny << " float" << endl;
   for(int j=0; j<ny; ++j)
      fout << grid_y[j] << " ";
   fout << endl;
   fout << "Z_COORDINATES " << nz << " float" << endl;
   for(int k=0; k<nz; ++k)
      fout << grid_z[k] << " ";
   fout << endl;

   fout << "POINT_DATA " << nx*ny*nz << endl;
   fout << "SCALARS density float" << endl;
   fout << "LOOKUP_TABLE default" << endl;
   // no need for k-loop since nk=1
   for(int i=0; i<nx; ++i)
      for(int j=0; j<ny; ++j)
      {
         for(int k=0; k<nz; ++k)
            fout << solution(i,j,k) << " ";
         fout << endl;
      }
   fout << "SCALARS density_exact float" << endl;
   fout << "LOOKUP_TABLE default" << endl;
   // no need for k-loop since nk=1
   for(int i=0; i<nx; ++i)
      for(int j=0; j<ny; ++j)
      {
         for(int k=0; k<nz; ++k)
            fout << solution_exact(i,j,k) << " ";
         fout << endl;
      }
   fout.close();
}
void vtk_anim_sol(vector<double> &grid_x,vector<double> &grid_y,
                  vector<double> &grid_z,
                  Array3D& solution, Array3D& solution_exact,
                   double t,
                  int time_step_number,
                  string filename)
{
  filename = filename+"_";
  filename = get_filename(filename,3,time_step_number);
  write_rectilinear_grid(grid_x, grid_y, grid_z, solution, solution_exact, t,
                         time_step_number, filename);
}