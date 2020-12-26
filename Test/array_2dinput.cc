#include "../include/array2d.h"

using namespace std;

int main()
{
  int ng = 1;
  double N = 2.0;
  Array2D solution(2,2, ng);
  solution = 0.0;
   //Will corner values in the second loop, i.e., with j
  //or we'd end up giving outdated values to corners.
  //i = -1, N
  solution(0,0) = 1.;
  int a[2];
  a[0] = 0, a[1] = 0;
  /*try
  {
    cout << solution(a);
  }*/
  //Thus, 2D input does NOT work
}