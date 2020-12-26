#include "../include/array2d.h"

using namespace std;

int main()
{
  int ng = 1;
  int N = 2;
  Array2D solution(2,2, ng);
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
  for (int k = -1; k<=N;k++)
  {
    solution(k,-1) = 0.0, solution(-1,k) = 0.0;
    solution(k,N) = 0.0 , solution(N,k) = 0.0;
  }
  solution.print_all(true);
}