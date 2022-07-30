#ifndef __ARRAY2D_H__
#define __ARRAY2D_H__

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip> //Used to define setw
//https://stdcxx.apache.org/doc/stdlibref/iomanip-h.html#:~:text=The%20header%20is%20part,the%20state%20of%20iostream%20objects.
#include <cassert>

using namespace std;

class Array2D
{
public:
   Array2D(); //Empty constructor declared
   Array2D(const int nx, const int ny,
           const int ng = 0/*Default value*/);
   void resize (const int nx, const int ny);
   void resize (const int nx, const int ny, const int ng);
   int a,b; //a,b are chosen such that
   //A(i,j) = u[a + i + j*b].
   //Actually, a = (2ng+nx+1)*ng, b = nx+2*ng
   //It is for optimizaiton that we are storing a,b separately.

   int sizex() const;
   int sizey() const;
   double  operator()(const int i, const int j) const;
   double& operator()(const int i, const int j);
   Array2D& operator= (const double scalar);
   Array2D& operator= (const Array2D& u);
   void update_fluff();
   //Overloads '<<', combining it with cout prints array without ghost cells
   //Question - Why is it inside the class?
   //'<<' overloaded as a friend so that it can access class variables.
   //The output is std::ostream format. When you call cout, it probably creates
   //and ostream object os which is fed to the arrows.
   //http://www.cplusplus.com/forum/beginner/99332/#msg534235
   //https://stackoverflow.com/questions/5508857/how-does-cout-actually-work
   //https://www.learncpp.com/cpp-tutorial/introduction-to-iostream-cout-cin-and-endl/
    friend std::ostream& operator<< (std::ostream&  os,
                                     const Array2D& A)
    {
      for(int i=0; i<A.sizex(); ++i)
      {
        for(int j=0; j<A.sizey(); ++j)
          os << std::setw(10) << A(i,j); //Sets width to be 10
        os << std::endl;
      }
      return os;
    }
    void print_all(bool label = false);

private:
  //We put ghost layers, which are extra columns/rows in our array
  //ng gives the number of ghost layers to be put on sides. If we put ng = 1,
  //we put one ghost layer each on left, right, top and bottom.
   int nx, ny, n, ng;
   std::vector<double> u;
};


#endif
