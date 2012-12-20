
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

//
// population.h -- functions, local variables, etc needed to simulate
// a single population
//

#ifndef _population_h_
#define _population_h_

#include <iostream>
#include <stddef.h>

class Individual;
class Generation;
class RNG;

class ResultBlock {
public:
  int ngen;			/* number of generations simulated */
  int nsur;			/* current size of the population */
};

class PopulationParamBlock {
public:
  double kmax;			/* K = mean carrying capacity */
  double sdk;			/*     std deviation of kmax */
  double rmax;			/* R = mean reproductive rate */
  double sdr;			/*     std deviation of rmax */
  double u;			/* mu = gametic mutation rate */
};

class Population {
public:
  Population(PopulationParamBlock *pb = NULL);
  ~Population();
  ResultBlock *run(int n);
  int step();
  friend std::ostream &operator <<(std::ostream &s, Population &p);
protected:
  Generation *cur;		/* current generation  */
  Generation *old;		/* past generation (only two age classes) */
  int ngen;			/* counts number of generations created */
  double maxkmax;		/* largest possible kmax */
  double kmax;			/* actual carrying capacity */
  double sdk;			/* std dev of kmax */
  double rmax;			/* reproductive rate */
  double sdr;			/* std dev of rmax */
  double u;			/* gametic mutation rate */
  RNG *rk;			/* RNG for kmax */
  RNG *rr;			/* RNG for rmax */
  RNG *ru;			/* RNG for new mutations */
  virtual int build_next_generation();
};

#endif  // _population_h_

