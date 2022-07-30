// Idea - (a,b,c) |-> a+nx*b+n*ny*c is a bijective map b/w
// (0,nx)*(0,ny)*(0,nz) and (0,nx*ny*nz)
#include <iostream>
#include "array3d.h"

using namespace std;

// Empty constructor
Array3D::Array3D ()
:
nx (0),
ny (0),
nz (0),
n  (0)
{
}

// Constructor based on size
Array3D::Array3D (const int nx, const int ny, const int nz)
:
nx (nx),
ny (ny),
nz (nz),
n  (nx*ny*nz),
u  (nx*ny*nz)
{
#if defined(DEBUG)
   if(nx < 0 || ny < 0 || nz < 0)
   {
      cout << "Dimensions cannot be negative" << endl;
      exit(0);
   }
#endif
}

// Change size of array
void Array3D::resize(const int nx1, const int ny1, const int nz1)
{
#if defined(DEBUG)
   if(nx1 < 0 || ny1 < 0)
   {
      cout << "Dimensions cannot be negative" << endl;
      exit(0);
   }
#endif
   nx = nx1;
   ny = ny1;
   nz = nz1;
   n  = nx1 * ny1 * nz1;
   u.resize (n);
}

// return number of rows, size of first index
int Array3D::sizex() const
{
   return nx;
}

// return number of columns, size of second index
int Array3D::sizey() const
{
   return ny;
}

int Array3D::sizez() const
{
   return nz;
}