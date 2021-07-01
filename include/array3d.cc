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

// Return value at (i,j), this is read only
double Array3D::operator() (const int i, const int j, const int k) const
{
#if defined(DEBUG)
   if(i < 0 || i > nx-1 || j < 0 || j > ny-1 || k < 0 || k > nz - 1)
   {
      cout << "Indices out of range" << endl;
      cout << "i, j = " << i << ", " << j << endl;
      exit(0);
   }
#endif
   return u[i + j*nx + k*nx*ny]; // TODO - Reorder so that k is the fastest index
}

// Return reference to (i,j), this can modify the value
double& Array3D::operator() (const int i, const int j, const int k)
{
#if defined(DEBUG)
   if(i < 0 || i > nx-1 || j < 0 || j > ny-1 || k < 0 || k > nz-1)
   {
      cout << "Indices out of range" << endl;
      cout << "i, j = " << i << ", " << j << endl;
      exit(0);
   }
#endif
   return u[i + j*nx + k*nx*ny];
}

// Set all elements to scalar value
Array3D& Array3D::operator= (const double scalar)
{
   for (unsigned int i=0; i<n; ++i)
      u[i] = scalar;
   return *this;
}

// Copy array a into this one
Array3D& Array3D::operator= (const Array3D& a)
{
#if defined(DEBUG)
   if(nx != a.sizex() || ny != a.sizey())
   {
      cout << "Array sizes do not match" << endl;
      exit(0);
   }
#endif
   u = a.u;
   return *this;
}