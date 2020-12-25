#include "array2d.h"

using namespace std;

int main()
{
  int ng = 1;
  double N = 2.0;
  Array2D solution(2,2, ng);
  solution = -1.0;
   //Will corner values in the second loop, i.e., with j
  //or we'd end up giving outdated values to corners.
  //i = -1, N
  solution(-1,-1) = 0.;
  solution(-1,0) = 1.;
  solution(-1,1) = 2.;
  solution(0,-1) = 3.;
  solution(1,-1) = 4.;
  for (int j = -1; j<N;j++) //not doing corners
  {
    solution(N-1,j) = solution(-1,j);
  //  solution(0,j)   = solution(N,j);
  }
  //j = -1, N
  for (int i = 0; i<N;i++) //doing corners
  {
    solution(i,N-1) = solution(i,-1);
   // solution(i,0)   = solution(i,N);
  }
  solution(1,1) +=1.0;
  cout << solution;
}