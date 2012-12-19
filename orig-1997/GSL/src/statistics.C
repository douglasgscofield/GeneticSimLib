
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

/*
  statistics.C -- parametric statistics calculation
*/

#include <stdlib.h>
#include <math.h>

#include "statistics.h"

const int NB=50;		// number of "bins" in the histogram
const int plotmin = 0;		// fix these (defined in Plot.h)
const int plotmax = NB;

Statistics::Statistics() {
  histogram = new int[NB];
  reset();
}

Statistics::~Statistics() {
  delete histogram;
}

void Statistics::reset() {
  int i;

  for (i=0; i<NB; i++)
    histogram[i] = 0;

  n = 0;
  sum = sumsq = 0.0;
  min = max = 0.0;
  mean = sd = cv = 0.0;
}

void Statistics::reset(double pmin, double pmax) {
  reset();
}

void Statistics::calculate() {
  double num = (double)n;
  double var;

  mean = sum/num;

  if (n > 1) {
    var = num * ((sumsq/num) - mean * mean) / (num - 1);
    sd = sqrt(var);
    cv = sd/mean;
  } 
  else
    sd = cv = 0.0;
}    

void Statistics::save(double x) {
  n += 1;
  if (n == 1)
    min = max = x;
  else {
    if (x > max)
      max = x;
    if (x < min)
      min = x;
  }
  sum += x;
  sumsq += (x * x);

  int i = floor((x-plotmin)*NB/(plotmax-plotmin));
  if ((i >= 0) && (i < NB))
    histogram[i] += 1;
}

StatsBlock *Statistics::get_results() {
  StatsBlock *res = new StatsBlock;

  calculate();

  res->n = n;
  res->mean = mean;
  res->sd = sd;
  res->cv = cv;
  res->min = min;
  res->max = max;

  return res;
}


