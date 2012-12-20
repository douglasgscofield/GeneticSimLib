
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

/*
  population.cpp -- implementation of class Population.  This class
  is basically just a place to collect all the odds 'n ends and
  run the simulation of a single population (the user interface
  coordinates simulation of multiple populations, possibly on
  different processors).
*/

#include <stddef.h>
#include <math.h>
#include <iostream>

#include "population.h"
#include "generation.h"
#include "individual.h"
#include "genome.h"
#include "rng.h"

// prototypes of local functions

Individual *mate(Individual *p1, Individual *p2, int nm);

// Constructor and destructor

Population::Population(PopulationParamBlock *pb) {
  kmax = pb->kmax;
  sdk = pb->sdk;
  rmax = pb->rmax;
  sdr = pb->sdr;
  u = pb->u * 2.0;			// convert gametic to zygotic rate

  maxkmax = kmax + 10.0*sdk;		// 10 sigma ought to be enough....
  ngen = 0;

  rk = new LogNormalRNG(kmax,sdk);
  ru = new PoissonRNG(u);
  rr = new LogNormalRNG(rmax,sdr);

  cur = NULL;
  old = new Generation;
  for (int i = 0; i < kmax; i++)	// should use sdk to figure first one?
    old->insert(new Individual);	// fix: add initial mutations
}

Population::~Population() {
  if (cur)
    delete cur;
  delete old;
  delete rk;
  delete ru;
  delete rr;
}

// run() -- keep building generations until one of them has no
// survivors (nsur = 0) or until the max number of steps has been reached
//
// Note: ngen is a state variable of this population object; the semantics
// of run() calls for us to continue the simulation using the current
// value of ngen....

ResultBlock *Population::run(int maxgen) {
  int nsur = old->size();
  
  while ((ngen < maxgen) && nsur) {
    nsur = build_next_generation();
    ngen += 1;
  }

  ResultBlock *r = new ResultBlock;
  r->ngen = ngen;
  r->nsur = nsur;
  return r;
}

// step() -- just do one step of the simulation, i.e. create one
// new generation.

int Population::step() {
  int nsur;
  int i;

  nsur = build_next_generation();
  ngen += 1;

  return nsur;
}

std::ostream &operator <<(std::ostream &s, Population &p) {
  s << "Population info: " << std::endl;
  s << "  cur = " << p.cur;
  if (p.cur) 
    s << " [n = " << p.cur->size() << "]";
  s << ", old = " << p.old;
  if (p.old)
    s << " [n = " << p.old->size() << "]";
  s << ",  ngen = " << p.ngen << ", kmax = " << p.kmax << ", maxkmax = " << p.maxkmax << std::endl;
  return s;
}

//    ------- local procedures --------

//
// build_next_generation() -- using rules for this particular kind
// of population, select parents from the current generation, mate
// them, and put surviving offspring in the new generation.
//
// This procedure implements monoecy:  adults are not differentiated
// by sex type, just grab two and combine them.
//
// This procedure is a candidate to move to a separate file or maybe
// a new virtual class type.
//

int Population::build_next_generation() {
  int  nparents;		/* size of old generation */
  int  nkids;			/* size of new generation */
  int  limit;			/* number of iterations in reproduction loop */
  int  target;			/* number of kids to create */
  Individual *kid;
  Individual *mom;
  Individual *dad;
  int  momx, dadx;
  fitness_t w;			// fitness of kid

  nparents = old->size();
  cur = new Generation;
  limit = nparents * rmax;
  target = floor(++(*rk) + 0.5);	// target is actual kmax for this generation
  if (target < 1)
    target = 1;
  nkids = 0;

  UniformRNG r(0.0,1.0);

  for (int i = 0; i < limit; i++) {
    momx = rword(nparents);
    dadx = rword(nparents);
    mom = (*old)[momx];
    dad = (*old)[dadx];
    kid = mate(mom,dad,++(*ru));		// create kid, add mutations
    if (kid->genes->fitness() > ++r) {
      cur->insert(kid);
      nkids += 1;
      if (nkids == target)
	break;
    }
    else {
      delete kid;
    }
  }

  delete old;			// update generations: toss the old one,
  old = cur;			// make the new one the current generation
  cur = NULL;

  return(nkids);
}

Individual *mate(Individual *p1, Individual *p2, int nm) {
  Individual *k = new Individual;

  k->genes->combine(p1->genes,p2->genes);
  k->genes->add_mutations(nm);
  return k;
}




