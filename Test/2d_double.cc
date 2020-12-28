#include<iostream>

using namespace std;

int main()
{
  double a[2];
  *a = (3.0,3.);
  if (a[0] == 0. && a[1] == 0.)
    cout <<"It works!\n";
  else 
  {
    cout <<"FAILURE!\n";
    cout <<"a[0] = "<<a[0] << ", a[1] = "<<a[1];
  }
}