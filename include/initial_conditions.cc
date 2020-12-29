#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <string>
#include <stdio.h>
#include <sys/time.h>

#include "/mnt/c/Users/arpit/Documents/GitHub/arpit_practise/include/array2d.h"
#include "/mnt/c/Users/arpit/Documents/GitHub/arpit_practise/include/vtk_anim.h"
using namespace std;

double interval_part(double x, double xmin, double xmax)
{
  if (x > xmax)
    return x - ceil((x - xmax) / (xmax - xmin)) * (xmax - xmin);
  else if (x < xmin)
    return x + ceil((xmin - x) / (xmax - xmin)) * (xmax - xmin);
  else
    return x;
}

double smooth_sine(double x, double xmin, double xmax)
{
  double value;
  value = sin(2.0 * M_PI * x / (xmax - xmin));
  return value;
}

//More testing needed before trying for different xmin,xmax
double hat_function(double x, double xmin, double xmax)
{
  double value;
  x = interval_part(x,xmin,xmax);
  if (x < xmin + (xmax - xmin) / 4.0 ||
      x > xmax - (xmax - xmin) / 4.0) //supported in (-1.5,1.5)
    value = 0.0;
  else if (x <= xmin+(xmax-xmin)/2.0)    //(x < 0.0)
    value = x -(xmin+ (xmax-xmin)/4.0);  //x + 0.5
  else if (x > xmax-(xmax-xmin)/2.0)     //x > 0.0
    value = (xmax-(xmax-xmin)/4.0) - x;  //0.5-x
  return value;
}

double step_function(double x, double xmin, double xmax)
{
  double value;
  x = interval_part(x,xmin,xmax);
  if (x < xmin + (xmax - xmin) / 4.0 || 
    x > xmax - (xmax - xmin) / 4.0)
    value = 0.0;
  else
    value = 1.0;
  return value;
}

double exp_func_25(double x, double xmin, double xmax)
{
  x = interval_part(x,xmin,xmax);
  return exp(-25*(x-0.25)*(x-0.25));
}

double exp_func_100(double x, double y, double xmin, double xmax, double ymin, double ymax)
{
  x = interval_part(x,xmin,xmax), y = interval_part(y,ymin,ymax);
  return 1.0 + exp(-100.0*((x-(xmin+0.25*(xmax-xmin)))*(x-(xmin+0.25*(xmax-xmin))) 
                            + (y-(0.5*(ymax+ymin)))*(y-(0.5*(ymax+ymin)))  ));
}

double cts_sine(double x, double xmin, double xmax)
{
    double value;
    x = interval_part(x,xmin,xmax);
    if (x < xmin + (xmax - xmin) / 4.0 ||
        x > xmax - (xmax - xmin) / 4.0)
        value = 0.0;
    else
        value = -sin(4.0 * M_PI * x / (xmax - xmin));
    return value;
}

class I_Functions
{
public:
  I_Functions();
  I_Functions(int initial_data_indicator,double xmin, double xmax);
  double value(double x, double y);

  //Gets exact solution, depending on advection speed
  double exact_value(double x, double y, double t, double u[2], 
                     bool constant = false);
  
  //Sets initial_data_indicator, xmin,xmax
  void set(int initial_data_indicator, double xmin, double xmax,
           double ymin, double ymax);
private:

  int initial_data_indicator;
  double xmin,xmax,ymin,ymax;
};

I_Functions::I_Functions(){}//Default constructor

I_Functions::I_Functions(int initial_data_indicator,
                           double xmin, double xmax):
                          initial_data_indicator(initial_data_indicator),
                          xmin(xmin),xmax(xmax),ymin(ymin),ymax(ymax)
{
  switch (initial_data_indicator)
  {
  case 0:
    cout <<"Smooth sine chosen for initial condition\n";
    break;
  case 1:
    cout <<"hat function chosen for initial condition\n";
    break;
  case 2:
    cout <<"discts step chosen for initial condition\n";
    break;
  case 3:
    cout <<"exp_25 chosen for initial condition\n";
    break;
  case 4:
    cout <<"exp_100 chosen for initial condition\n";
    break;
  case 5:
    cout <<"cts_sine chosen for initial condition \n";
    break;
  default:
    cout << "You entered the wrong initial_data_indicator ";
    assert(false);
  }
}

void I_Functions::set(int initial_data_indicator0, double xmin0, double xmax0,
                      double ymin0, double ymax0)
{
  initial_data_indicator = initial_data_indicator0;
  if (initial_data_indicator0 == 0)
    cout << "WARNING - smooth_sine doesn't work for variable coefficients\n";
  xmin = xmin0, xmax = xmax0;
  ymin = ymin0, ymax = ymax0;
}

//BUG - This way only works for (xmin,xmax) = (ymin,ymax).
double I_Functions::value(double x, double y)
{
  double val, radius;
  radius = (x-0.5)*(x-0.5)*y*y;//centre of function, used in bump function.
  switch (initial_data_indicator)
  {
  case 0:
    return sin(2.0 * M_PI * x / (xmax - xmin))* sin(2.0 * M_PI * y / (ymax - ymin));
    break;
  case 1:
    return hat_function(x,xmin,xmax)*hat_function(y,ymin,ymax);
    break;
  case 2:
    return step_function(x,xmin,xmax)*step_function(y,ymin,ymax);
    break;
  case 3:
    return exp_func_25(x,xmin,xmax)*exp_func_25(y,ymin,ymax);
    break;
  case 4:
    return exp_func_100(x,y,xmin,xmax,ymin,ymax);
    break;
  case 5:
    return cts_sine(x,xmin,xmax)*cts_sine(y,ymin,ymax);
    break;
  default:
    cout << "You entered the wrong initial_data_indicator ";
    assert(false);
  }
}

double I_Functions::exact_value(double x, double y,double t, double u[2], bool constant)
{
  if (constant == true)
  {
    return value(x-u[0]*t,y-u[1]*t);
  }
  else
  {
    return value(x*cos(t)+y*sin(t),-x*sin(t)+y*cos(t));
  }
}


