#ifndef __VTK_ANIM_H__
#define __VTK_ANIM_H__
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <cassert>

#include "array2d.h"

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
                            Array2D &solution,
                            double t,
                            int c, //Cycle number
                            string filename)
{
   const int nx = solution.sizex(),ny = solution.sizey();
   if ( (grid_x.size() != nx) || (grid_y.size() != ny))
    {
      cout << "Grid and solution sizes don't match." <<endl;
      cout << "grid.sizex = "<< grid_x.size() << ", solution.sizex = "<<nx;
      cout << "grid.sizey = "<< grid_y.size() << ", solution.sizey = "<<ny;
      assert(false);
    }
   int nz = 1; // We have a 2d grid
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
   fout << 0.0 << endl;

   fout << "POINT_DATA " << nx*ny*nz << endl;
   fout << "SCALARS density float" << endl;
   fout << "LOOKUP_TABLE default" << endl;
   // no need for k-loop since nk=1
   for(int j=0; j<ny; ++j)
   {
      for(int i=0; i<nx; ++i)
         fout << solution(i,j) << " ";
      fout << endl;
   }
   /*
   fout << "SCALARS density_exact float" << endl;
   fout << "LOOKUP_TABLE default" << endl;
   // no need for k-loop since nk=1
   for(int j=0; j<ny; ++j)
   {
      for(int i=0; i<nx; ++i)
         fout << solution_exact(i,j) << " ";
      fout << endl;
   }
   fout.close();
   */

   //cout << "Wrote Cartesian grid into " << filename << endl;
}
void vtk_anim_sol(vector<double> &grid_x,vector<double> &grid_y,
                  Array2D& solution,
                   double t,
                  double time_step_number,
                  string filename)
{
  //cout << "Saving rectilinear grid\n";
  /*
  const int nx = solution.sizex(), ny = solution.sizey();
  const double xmin = -1.0, xmax = 1.0;
  const double ymin = -1.0, ymax = 1.0;
  const double dx = (xmax-xmin)/(nx);
  const double dy = (ymax-ymin)/(ny);
  */
  /*double *x, *y;

  x = new double[nx];
  y = new double[ny];

  for(int i=0; i<nx; ++i)
    x[i] = xmin + i*dx;

  for(int j=0; j<ny; ++j)
    y[j] = ymin + j*dy;
  */
  /*
  var = new double*[nx];
  for(int i=0; i<nx; ++i)
    var[i] = new double[ny];
  */
  //double u = 1.0, v = 1.0;
  //double t = 0.0, dt = 0.1, Tf = 10.0;
  /*
  for(int i=0; i<nx; ++i)
      for(int j=0; j<ny; ++j)
      {
        //double xx = x[i] - u*t;
        //double yy = y[j] - v*t;
        var[i][j] = solution(i,j);
      }
  */
  filename = filename+"_";
  filename = get_filename(filename,3,time_step_number);
  write_rectilinear_grid(grid_x, grid_y, solution, t, time_step_number, filename);
  //t += dt;
  //++c;

   // Deallocate memory here
   /*delete[] x;
   delete[] y;*/
   /*for(int i=0; i<nx; ++i)
      delete[] var[i];
   delete[] var;*/

}
#endif
