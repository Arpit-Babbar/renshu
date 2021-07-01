#include "../../include/array3d.h"
#include "../../include/array3d.cc"
// Test whether Array3D can be given as a function argument
using namespace std;
void tester(Array3D solution[])
{
  cout << solution[0]<<endl;
  cout << solution[1]<<endl;
}
int main()
{
  int nx = 2, ny = 2, nz = 3;
  int levels = 2;
  Array3D solution[levels]; // solution old and solution new
  solution[0].resize(nx,ny,nz), solution[1].resize(nx,ny,nz);
  double iter = 0.0;
  for(int i=0;i<nx;i++)
    for(int j=0;j<ny;j++)
      for(int k=0;k<nz;k++)
      {
        solution[0](i,j,k) = iter;
        solution[1](i,j,k) = 10+iter;
        iter++;
      }
  tester(solution);

}