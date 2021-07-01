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
   double  operator()(const int i, const int j, const int k) const; // read value
   double& operator()(const int i, const int j, const int k); // read/write value
   Array3D& operator= (const double scalar); // set array to constant
   Array3D& operator= (const Array3D& u);    // set array to another array
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
