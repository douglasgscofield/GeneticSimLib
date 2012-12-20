
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

//
// Derived classes for Coho Salmon simulation
//
// The classes defined in this file are derived from classes in the 
// genetics simulation library:
//
//    Class		Base			Description
//
//    Coho		Individual		A single coho salmon
//    CohoPopulation	Population		Three generations of salmon
//    CohoParamBlock	PopulationParamBlock	Extra parameters needed by
//						a coho population

#include <stddef.h>
#include "gsl.hpp"

// The extra data member that defines a salmon object is a "return_year"
// that determines when the individual will return to spawn.  Both the sex 
// and return year of a new object are set at birth and never changed.

class Coho : public Individual {
public:
  Coho(double pf, double pm, double pj);
  int return_year();
protected:
  int d_return;
};

// The additional parameters in the coho salmon simulation are the 
// sex ratio and return year.  The classes are females, which always
// return to spawn in their third year; males that return as three-year
// olds, and jacks, which are males which return as two-year olds.
// The parameters are the probability that a new individual falls
// into that class, and are used by the Coho constructor to assign
// the sex and return year of the new individual.

class CohoParamBlock : public PopulationParamBlock {
public:
  double pfemale;		/* pct of offspring that are females */
  double pmale;			/* pct that are males returning in 3 years */
  double pjack;			/* pct that are males returning in 2 years */
};

// A CohoPopulation, derived from Population, inherits the state variables, 
// initialization procedures, and other routines.  Instead of just two
// generations (cur and old) it keeps four:  three generations of adults and
// a single new generation.  This class overrides the definition of
// "build_next_generation" to use two previous generations to produce
// the next one.

static const int NYEARS = 4;	// number of generations to manage

class CohoPopulation : public Population {
public:
  CohoPopulation(CohoParamBlock *pb, int *as, char *id);
  ~CohoPopulation();
  friend std::ostream &operator <<(std::ostream &s, CohoPopulation &p);
private:
  int build_next_generation();
  void new_generation();
  void build_parent_generation();
  void build_mating_groups();
  void create_offspring();
  Generation *m[NYEARS];	// m[i] is set of i-year old males
  Generation *f[NYEARS];	// f[i] is set of i-year old females
  Generation *j[NYEARS];	// j[i] is set of i-year old jacks
  ISet *mp;			// male parents of next generation
  ISet *fp;			// female parents of next generation
  double pfemale;
  double pmale;
  double pjack;
  // items saved for printing at end of each generation
  char *pid;			// population id
  int genx;			// generation number
  int px;			// number of returns
  int pmx;			// number of male parents
  int pfx;			// number of female parents
  float rx;			// reproductive rate used to build this generation
  int lx;			// limit (number of pairings)
  int kx;			// carrying capacity of this generation (target)
  int mx;			// number of males, this generation
  int jx;			// number of jacks, this generation
  int fx;			// number of females, this generation
  fitness_t wx;			// mean fitness, survivors
  fitness_t vwx;		// std dev fitness, survivors				     
  fitness_t awx; 		// mean fitness, all progeny
};

typedef enum {
  COHO_INIT_ERROR,
  COHO_NO_COHO,
  COHO_NO_NEW_GENERATION,
  COHO_NO_PARENTS,
  COHO_NO_CHILD
} COHO_ALLOC_ERROR;
