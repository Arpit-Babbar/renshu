#ifndef __INITIAL_CONDITION_H__
#define __INITIAL_CONDITION_H__

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <cassert>
#include <string>
#include <stdio.h>
#include <sys/time.h>
using namespace std;

double interval_part(double x, double xmin, double xmax);
double smooth_sine(double x, double xmin, double xmax);
//More testing needed before trying for different xmin,xmax
double hat_function(double x, double xmin, double xmax);
double step_function(double x, double xmin, double xmax);

double exp_func_25(double x, double xmin, double xmax);
double exp_func_50(double x, double y, double xmin, double xmax, double ymin, double ymax);

double cts_sine(double x, double xmin, double xmax);

class I_Functions
{
public:
  I_Functions();
  I_Functions(int initial_data_indicator,double xmin, double xmax,
              double ymin, double ymax);
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
#endif
