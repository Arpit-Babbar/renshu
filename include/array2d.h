#ifndef __ARRAY2D_H__
#define __ARRAY2D_H__

#include <vector>

using namespace std;

class Array2D
{
public:
   Array2D();
   Array2D(const unsigned int nx, const unsigned int ny/*,
           const unsigned int ng = 0*//*Default value*/); 
   void resize (const unsigned int nx, const unsigned int ny);
   unsigned int sizex() const;
   unsigned int sizey() const;
   double  operator()(const unsigned int i, const unsigned int j) const;
   double& operator()(const unsigned int i, const unsigned int j);
   Array2D& operator= (const double scalar);
   Array2D& operator= (const Array2D& u);
   
private:
   /*const */unsigned int nx, ny, n/*, ng*/;
   std::vector<double> u;
};


// Empty constructor
Array2D::Array2D ()
:
nx (0),
ny (0),
/*ng (0),*/
n  (0)
{
}
//We can specify a layer of ghost cells 
//Adding a print function.
// Constructor based on size
Array2D::Array2D (const unsigned int nx, const unsigned int ny/*,
                  const unsigned int ng*/)
:
nx (nx),
ny (ny),
/*ng (ng),*/
n  ((nx/*+2*ng*/)*(ny/*+2*ng*/)),
u  ((nx/*+2*ng*/)*(ny/*+2*ng*/))
{
}

//Overload constructor.



// Change size of array
void Array2D::resize(const unsigned int nx1, const unsigned int ny1)
{
   nx = nx1;
   ny = ny1;
   n  = nx1 * ny1;
   u.resize (n);
}

// return number of rows, size of first index
unsigned int Array2D::sizex() const
{
   return nx;
}

// return number of columns, size of second index
unsigned int Array2D::sizey() const
{
   return ny;
}

// Return value at (i,j), this is read only
double Array2D::operator() (const unsigned int i, const unsigned int j) const
{
  if ((i>=nx) || (j>=ny))
   {
   cout << "Attempt to access non-existent array entries"<<endl;
   cout << "Array has rows,columns of sizes " <<  nx <<","<< ny << endl;
   cout << "Tried to access row, column position " << i <<"," <<j << endl;
   assert(false);
   }
   return u[i + j*nx];
}

// Return reference to (i,j), this can modify the value
double& Array2D::operator() (const unsigned int i, const unsigned int j)
{
#ifdef DEBUG /* g++ -o output main.cc -DDEBUG*/
   if ((i>=nx) || (j>=ny)) //Move to debug. Debugging option in makefile.
   {
   cout << "Attempt to access non-existent array entries"<<endl;
   cout << "Array has rows,columns of sizes " <<  nx <<","<< ny << endl;
   cout << "Tried to access row, column position " << i <<"," <<j << endl;
   assert(false);
   }
#endif
   return u[i + j*nx];
}

// Set all elements to scalar value
Array2D& Array2D::operator= (const double scalar)
{
   for (unsigned int i=0; i<n; ++i)
      u[i] = scalar;
   return *this;
}

// Copy array a into this one
Array2D& Array2D::operator= (const Array2D& a)
{
   u = a.u;
   return *this;
}


#endif
