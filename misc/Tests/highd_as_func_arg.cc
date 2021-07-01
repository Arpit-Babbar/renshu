#include <iostream>
#include <cmath>
#define p_dim 3

using namespace std;

void print_array(int arr[2][p_dim])
{
   cout<<"left limit values"<<endl;
   for(int i = 0; i<p_dim; i++)
      cout<<arr[0][i] << endl;
   cout<<"right limit values"<<endl;
   for(int i = 0; i<p_dim; i++)
      cout<<arr[1][i] << endl;
}

// Test whether Array3D can be given as a function argument
using namespace std;

int main()
{
   int udim[2][p_dim]={0};
   print_array(udim);
}