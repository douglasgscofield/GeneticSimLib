
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

/*
  coho.C -- implementation of classes used in the Coho salmon genetic
  simulations.  Classes defined here are derived from the Individual, 
  Population, and PopulationParamBlock classes of the Genetic
  Simulation Library.  Each derived class inherits all the
  state variables and public functions of the parent class and
  extends it with characteristics of coho salmon -- e.g. the
  Population class implements an age structure such that parents
  used to create a new generation are drawn from 3-year old females
  and 2- and 3-year old males according to proportions specified
  in a CohoParamBlock object.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <math.h>
#include <iostream.h>
#include <iomanip.h>
#include <time.h>
#include <unistd.h>
#include <sys/mman.h>

#include "gsl.h"
#include "coho.h"

void delete_set(ISet **p);

// Class Coho is derived from Individual.  The constructor for a Coho
// object assigns the sex and return year at birth.  The simulation parameters 
// are 'pfemale', 'pmale', and 'pjack', which divide the population into three 
// groups.  All members of 'male' and 'jack' are males, the rest are females.  
// All members of 'jack' return at 2 years, the rest return in 3 years.

Coho::Coho(double pf, double pm, double pj) {
  double x = runiform(0.0,1.0);

  d_sex = (x < pf) ? FEMALE : MALE;
  d_return = (x < (pf + pm)) ? 3 : 2;    
}


int Coho::return_year() {
  return d_return;
}

// Every instance of a CohoPopulation will have 4 generations,
// named g[0] to g[3] (to keep it general, the .h file defines
// a constant NYEARS = 4).  g[0] will be the new generation built by
// a call to build_next_generation, and g[n], 1 <= n <= 3, will be the
// generation of n-year olds.  

CohoPopulation::CohoPopulation(CohoParamBlock *pb, int *as, char *id) : Population(pb) {
  int i, k;
  UniformRNG r(0.0,1.0);
  Coho *x;
  double r0;

  m[0] = f[0] = j[0] = NULL;
  sdr = pb->sdr;

  pfemale = pb->pfemale;
  pmale = pb->pmale;
  pjack = pb->pjack;

  // The parent class (Population) creates LogNormal random number generators
  // for R and K.  In this simulation, we want to use LogNormalLog RNGs (since
  // the mean and std dev of the underlying Gaussians are input parameters).
  // Also, we scale R by exp(u) to account for initial mutations.

  delete rk;
  rk = new LogNormalLogRNG(kmax,sdk);

  delete rr;
  r0 = rmax + u;			// mean draw will be exp(rmax)*exp(u)
  rr = new LogNormalLogRNG(r0,sdr);

  for (i = 1; i < NYEARS; i++) {
    m[i] = new Generation;
    f[i] = new Generation;
    j[i] = new Generation;
    for (k = 0; k < as[i-1]; k++) {
      if (x = new Coho(pfemale,pmale,pjack)) {
	if (x->get_sex() == MALE)
	  if (x->return_year() == 3)
	    m[i]->insert(x);
	  else
	    j[i]->insert(x);
	else
	  f[i]->insert(x);
	x->genes->combine(NULL,NULL); 	// assigns random fitness
      }
      else
	throw COHO_NO_COHO;
    }
  }

  mp = fp = NULL;

  // Toss out any 3-year old jacks (all jacks would have returned
  // the previous year)

  delete j[3];

  // Create the population id and initialize other params
  // printed at the end of each generation:

  static int cpcount = 1;

  pid = new char[20];
  sprintf(pid,"%s:%d",id,cpcount++);

  genx = 0;
  px = pmx = pfx = lx = kx = mx = fx = jx = 0;
  awx = vwx = wx = rx = 0.0;
}

CohoPopulation::~CohoPopulation() {
  int i;

  for (i = 0; i < NYEARS; i++) {
    delete m[i];
    delete f[i];
  }

  for (i = 0; i < NYEARS-1; i++)
    delete j[i];

  if (mp)
    delete_set(&mp);
  if (fp)
    delete_set(&fp);

  delete pid;
}

ostream &operator <<(ostream &s, CohoPopulation &p) {
  int i, k;
  Generation *g;
  ISet *ip;

  for (i = 0; i < NYEARS; i++) {
    s << i << "-year olds:" << endl;
    s << "  M: ";
    if (g = p.m[i]) {
      for (k = 0; k < g->size(); k++)
	s << setw(4) << (*g)[k]->get_id();
    }
    s << endl;
    s << "  F: ";
    if (g = p.f[i]) {
      for (k = 0; k < g->size(); k++)
	s << setw(4) << (*g)[k]->get_id();
    }
    s << endl;
    s << "  J: ";
    if (g = p.j[i]) {
      for (k = 0; k < g->size(); k++)
	s << setw(4) << (*g)[k]->get_id();
    }
    s << endl;
  }

  s << endl;
  if (p.mp && p.fp) {
    s << "Parent set:" << endl;
    s << "  M: ";
    ip = p.mp;
    for (k = 0; k < ip->size(); k++)
      s << setw(4) << (*ip)[k]->get_id();
    s << endl;
    s << "  F: ";
    ip = p.fp;
    for (k = 0; k < ip->size(); k++)
      s << setw(4) << (*ip)[k]->get_id();
    s << endl;
  }
  s << endl;

  return s;
}


// build_next_generation()
//
// This procedure simulates the actions in a complete year.  It creates
// a new generation of offspring, adds them to the overall population,
// and then ages all individuals by one year.
//
// Start a new generation by first selecting the breeding population.  This
// is done by "sampling without replacement" -- move a number of individuals
// from previous generations into the parent sets fp and mp.  Then pair them
// up and create offspring, putting survivors into m[0] and f[0].  NOTE: if
// one parent set is empty (e.g. no females selected for mating) then m[0]
// and f[0] will both be empty sets.

int CohoPopulation::build_next_generation() {
  int i;
  int n;

  time_t start_time = clock();

  new_generation();			// clear out parent sets, choose params, etc
  build_parent_generation();		// puts parent groups in 'mp', 'fp'
  create_offspring();			// do matings and put survivors in m[0] and f[0]

  delete m[NYEARS-1];			// done with year 3 males and females...
  delete f[NYEARS-1];
  delete j[NYEARS-2];			// and year 2 jacks

  for (i = NYEARS-2; i >= 0; i--) { 	// age the remaining generations by one year
    m[i+1] = m[i];	
    f[i+1] = f[i];
    j[i+1] = j[i];
  }

  m[0] = f[0] = j[0] = NULL;
  
  time_t end_time = clock();
  double t = ((double)end_time - (double)start_time)/CLOCKS_PER_SEC;

  // Output complete information about this generation:

  genx += 1;
  const int obufsize = 100;
  float cvx = wx ? vwx/wx : 0.0;
  static char obuf[obufsize];
  sprintf(obuf,"%s G%03d %5d %5d %5d %5d %6.2f %5d %5d %5d %5d %5.3f %5.3f %5.3f %5.3f %5.3f",
	  pid, genx, px, kx, pmx, pfx, rx, lx, fx, mx, jx, awx, wx, vwx, cvx, t);
  cout << obuf << endl;

  // The simulation is now "between generations" -- everything is
  // set for the next call to build_next_generations:  all generations
  // have aged one year, the most recent generation is in m[1], j[1], and 
  // f[1], and m[0], j[0], and f[0] are ready to receive the next generation.
  // Return the total number of individuals in the population.

  n = 0;
  for (i = 1; i < NYEARS; i++)
    n += m[i]->size() + f[i]->size() + j[i]->size();

  return n;
}

// new_generation() -- a local procedure that re-initializes the object
// state variables at the beginning of each new generation.  We could
// deallocate these things when we are done with them at the end of
// build_next_generation, but they are left in state variables in case
// a demo program wants to see them (via a call to the << operator)
// after one generation is built and before the next one is started.
//
// The other initialization is to choose the size of the parent set
// and their reproductive rates.

void CohoPopulation::new_generation() {
  if (mp) {
    delete_set(&mp);
    mp = NULL;
  }
  if (fp) {
    delete_set(&fp);
    fp = NULL;
  }

  if ((m[0] = new Generation) && (f[0] = new Generation) && (j[0] = new Generation))
    return;

  throw COHO_NO_NEW_GENERATION;
}


// Local procedure for this type of population -- move individuals
// out of previous generations and build a set of parents that will
// be used to produce the new generation.
//
// The first step is to create a set of "returns" consisting of the
// entire set of 3-year old females, jacks from the set of 2-year old
// males, and the remaining 3-year old males.  Then use a log-normal
// random number generator to set the size of the parent set, and
// finally choose that many individuals at random from the pool.
// Male parents go to the set mp, female parents to fp.

void CohoPopulation::build_parent_generation() {
  Coho *ip;
  int target;
  int i, k;
  int np, nmp, njp, nfp;
  Generation *from;
  ISet *to;

  fp = new ISet;
  mp = new ISet;

  if ((fp == NULL) || (mp == NULL))
    throw COHO_NO_PARENTS;

  nfp = f[3]->size();
  nmp = m[3]->size();
  njp = j[2]->size();
  px = np = nfp + nmp + njp;
  kx = target = int(++(*rk)+0.5);

  if (np < target)
    target = np;

  while (target) {
    k = rword(np);
    if (k < nfp) {
      from = f[3];
      to = fp;
    }
    else if (k < (nfp+nmp)) {
      from = m[3];
      to = mp;
    }
    else {
      from = j[2];
      to = mp;
    }
    if (from->size()) {
      to->insert(from->remove(0));
      target -= 1;
    }
  }
}


// Create the new generation using the parent sets 'mp' and 'fp'.  Use 
// random selection of mates.

static int gencount = 0;

void CohoPopulation::create_offspring() {
  int limit;
  int i, k, nm, nf;
  Coho *child;
  UniformRNG r(0.0,1.0);
  float ri;
  fitness_t w, mw, vw, amw;

  pmx = nm = mp->size();
  pfx = nf = fp->size();

  rx = ri = ++(*rr);		// get sample from log normal
  ri = ri / pfemale;		// scale by percentage of females
  limit = ri * nf;		// this is the number of offspring to make
  lx = limit;

  amw = mw = vw = 0.0;
  mx = jx = fx = 0;

  //  static int gcount = 1;			// @@@
  //  cout << "Gen " << gcount++ << endl; 		// @@@

  //  mlockall(MCL_CURRENT | MCL_FUTURE);

  if ((nm > 0) && (nf > 0)) {
    for (i = 0; i < limit; i++) {
      //      p1 = (Coho *)(*fp)[rword(nf)];
      //      p2 = (Coho *)(*mp)[rword(nm)];
      if ((child = new Coho(pfemale,pmale,pjack)) == NULL)
	throw COHO_NO_CHILD;
      //      child->genes->combine(p1->genes,p2->genes);
      w = child->genes->combine(NULL,NULL);
      //      w = child->genes->add_mutations(++(*ru));
      amw += w;
      if (w > ++r) {
	// cout << w << endl;			// @@@
	mw += w;
	vw += w*w;
	if (child->get_sex() == MALE)
	  if (child->return_year() == 2) {
	    j[0]->insert(child);
	    jx += 1;
	  }
	  else {
	    m[0]->insert(child);
	    mx += 1;
	  }
	else {
	  f[0]->insert(child);
	  fx += 1;
	}
      }
      else {
	delete child;
      }
    }
  }

  //  munlockall();

  int n = mx+fx+jx;

  if (n > 1) {
    wx = mw/n;
    vwx = sqrt(n*(vw/n - wx*wx)/(n-1));
  }
  else {
    wx = 0.0;
    vwx = 0.0;
  }

  if (limit)
    awx = amw/limit;
  else
    awx = 0.0;  
}

// Delete all the individuals in an ISet and then the set itself

void delete_set(ISet **p) {
  ISet *ps = *p;

  while (ps->size())
    delete ps->remove(0);
  delete ps;
  *p = NULL;
}

