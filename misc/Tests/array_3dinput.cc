#include "../../include/array3d.h"
#include "../../include/array3d.cc"

using namespace std;

int main()
{
  double N = 2.0;
  int nx = 2, ny = 2, nz = 3;
  Array3D solution(nx,ny,nz);
  solution = 0.0;
  double iter = 0.0;
  for(int i=0;i<nx;i++)
    for(int j=0;j<ny;j++)
      for(int k=0;k<nz;k++)
      {
        printf("i,j,k = %d,%d,%d\n", i,j,k);
        solution(i,j,k) = iter;
        iter++;
      }
  cout << solution;
}