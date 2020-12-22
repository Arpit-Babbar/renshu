#include "array2d.h"

using namespace std;

int main()
{
  Array2D test(3.0,4.0);
  test = 0.0;
  double c = 0.0;
  for (int i = 0; i < test.sizex(); i++)
    for (int j = 0; j<test.sizey(); j++)
    {
      test(i,j) = c;
      c = c + 1.0;
    }
  cout << test << endl;
}