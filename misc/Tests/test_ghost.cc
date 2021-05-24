#include "../fv2d_var_coeff/array2d.h"
#include "../fv2d_var_coeff/array2d.cc"

using namespace std;

int main()
{
  int ng = 1;
  int N = 2;
  Array2D solution(2,3, ng);
  solution = 0.0;
   //Will corner values in the second loop, i.e., with j
  //or we'd end up giving outdated values to corners.
  //i = -1, N
  solution(-1,-1) = 1.;
  solution(-1,0) = 2.;
  solution(0,-1) = 3.;
  solution(1,-1) = 4.;
  solution(-1,1) = 5.0;
  solution(2,2)=6.0;
  solution(2,1) = 7.0;
  solution(1,2) = 8.0;
  solution.print_all(true);
  
  solution.update_fluff();
  cout <<"Updated matrix is "<<endl;

  solution.print_all(true);
  
  solution = 0.0;
  double c = 0;
  cout << "Here's another empty matrix"<<endl;
  for (int i = 0; i<solution.sizex();i++)
    for (int j = 0; j<solution.sizey();j++)
    {
      solution(i,j) = c;
      c++;
    }
  solution.print_all(true);  
  cout <<"After update_fluff, updated matrix is "<<endl;
  solution.update_fluff();
  solution.print_all(true);
}