
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

//
// statistics.h -- calculate parametric statistics 
//

#ifndef _statistics_h_
#define _statistics_h_

class ResultBlock;

class StatsBlock {
public:
  int n;			/* N */
  double mean;			/* M */
  double sd;			/* standard deviation */
  double cv;			/* coefficient of variation */
  int min;
  int max;
};

class Statistics {
public:
  Statistics();
  ~Statistics();
  void reset();
  void reset(double min, double maxm);
  void save(double x);
  StatsBlock *get_results();
protected:
  void calculate();
  int *histogram;
  int n;
  double sum;
  double sumsq;
  double max;
  double min;
  double mean;
  double sd;
  double cv;
};

#endif  // _statistics_h_

