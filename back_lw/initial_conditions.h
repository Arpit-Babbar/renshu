#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include <cassert>
using namespace std;

class Initial_Data
{
  public:
  void set_initial_data(double x_min0, double x_max0, string initial_data_indicator0);
  double value(double x);
  private:
  double x_min,x_max;
  string initial_data_indicator;
  double interval_part(double x);
  double hat_function(double grid_point);
  double step_function(double grid_point);
  double sine_wave(double grid_point);
};


