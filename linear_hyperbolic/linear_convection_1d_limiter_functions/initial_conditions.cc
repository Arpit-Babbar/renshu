#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include <cassert>
#include "initial_conditions.h"
using namespace std;

double Initial_Data::value(double x)
{
  if (initial_data_indicator=="hat")
  {
    double result = hat_function(x);
    return result;
  }
  else if (initial_data_indicator=="smooth_sine")
  {
    double result = sin(2.0 * M_PI * x / (x_max - x_min));
    return result;
  }
  else if (initial_data_indicator=="step")
  {
    double result = step_function(x);
    return result;
  }
  else if (initial_data_indicator=="cts_sine")
  {
    double result = sine_wave(x);
    return result;
  }
  else if (initial_data_indicator=="exp_func_100")
  {
    double result = exp_func_100(x);
    return result;
  }
  else 
  {
    cout <<"Incorrect initial_data_indicator"<<endl;
    cout <<"You put "<<initial_data_indicator<<endl;
    cout <<"Your choices are "<<endl;
    cout << "hat, smooth_sine, step, cts_sine, exp_func_100 ." <<endl;
    assert(false);
  }
}
void Initial_Data::set_initial_data(double x_min0, double x_max0, string initial_data_indicator0)
{
  x_min= x_min0;
  x_max = x_max0;
  initial_data_indicator = initial_data_indicator0;
}
double Initial_Data::interval_part(double x)
{
  if (x > x_max)
      return x - ceil((x - x_max) / (x_max - x_min)) * (x_max - x_min);
  else if (x < x_min)
      return x + ceil((x_min - x) / (x_max - x_min)) * (x_max - x_min);
  else
      return x;
}
double Initial_Data::hat_function(double grid_point)
{
  double result;
  grid_point = interval_part(grid_point);
  if (grid_point <= -0.5 || grid_point >= 0.5)
    result = 0.0;
  else 
  {
    if (grid_point <= 0.0)
      result = grid_point + 0.5; //x+1/2
    if (grid_point >  0.0)
      result = 0.5 - grid_point; //x-1/2
  }
  return result;
}

double Initial_Data::step_function(double grid_point)
{
  double value;
  grid_point = interval_part(grid_point);
  if (grid_point < x_min + (x_max - x_min) / 4.0 ||
    grid_point > x_max - (x_max - x_min) / 4.0)
    value = 0.0;
  else
    value = 1.0;
  return value;
}

double Initial_Data::sine_wave(double grid_point)
{
  double value;
  grid_point = interval_part(grid_point);
  if (grid_point < x_min + (x_max - x_min) / 4.0 ||
    grid_point > x_max - (x_max - x_min) / 4.0)
    value = 0.0;
  else
    value = -sin(4.0 * M_PI * grid_point / (x_max - x_min));
  return value;
}

double Initial_Data::exp_func_100(double grid_point)
{
  double value;
  grid_point = interval_part(grid_point);
  value = exp(-100.0*grid_point*grid_point);
  return value;
}

