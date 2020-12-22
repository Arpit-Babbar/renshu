#include "array2d.h"

using namespace std;

int main()
{
  Array2D test(1.0,2.0,1.0);
  test = 0.0;
  double c = 0.0;
  for (int i = -1; i <= test.sizex(); i++)
    for (int j = -1; j<=test.sizey(); j++)
    {
      test(i,j) = c;
      c = c + 1.0;
    }
  cout << test << endl;
}