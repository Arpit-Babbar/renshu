#ifndef __ARRAY2D_H__
#define __ARRAY2D_H__

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip> //Used to define setw
//https://stdcxx.apache.org/doc/stdlibref/iomanip-h.html#:~:text=The%20header%20is%20part,the%20state%20of%20iostream%20objects.
#include <cassert>

using namespace std;

// TODO - Add ghost cells

class Array3D
{
public:
   Array3D(); //Empty constructor declared
   Array3D(const int nx, const int ny, const int nz);
   void resize (const int nx, const int ny, const int nz);
   int sizex() const;
   int sizey() const;
   int sizez() const;
   // Return value at (i,j), this is read only
  double operator() (const int i, const int j, const int k) const
  {
  #if defined(DEBUG)
    if(i < 0 || i > nx-1 || j < 0 || j > ny-1 || k < 0 || k > nz - 1)
    {
        cout << "Indices out of range" << endl;
        cout << "i, j = " << i << ", " << j << endl;
        exit(0);
    }
  #endif
    return u[i*ny*nz + j*nz + k]; // TODO - Reorder so that k is the fastest index
  }
   // Return reference to (i,j), this can modify the value
  double& operator() (const int i, const int j, const int k)
  {
  #if defined(DEBUG)
    if(i < 0 || i > nx-1 || j < 0 || j > ny-1 || k < 0 || k > nz-1)
    {
        cout << "Indices out of range" << endl;
        cout << "i, j = " << i << ", " << j << endl;
        exit(0);
    }
  #endif
    return u[i*ny*nz + j*nz + k];
  }
   // Set all elements to scalar value
  Array3D& operator= (const double scalar)
  {
    for (int i=0; i<n; ++i)
        u[i] = scalar;
    return *this;
  }
   // Copy array a into this one
  Array3D& operator= (const Array3D& a)
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
    friend std::ostream& operator<< (std::ostream&  os,
                                     const Array3D& A)
    {
      for(int i=0; i<A.sizex(); ++i)
      {
         for(int j=0; j<A.sizey(); ++j)
           for (int k=0; k<A.sizez(); ++k)
             os << std::setw(10) << A(i,j,k); //Sets width to be 10
        os << std::endl;
      }
      return os;
    }
private:
   int nx, ny, nz, n;
   std::vector<double> u;
};


#endif
