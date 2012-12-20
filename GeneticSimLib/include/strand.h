
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

//
// strand.h -- used by SparseGenome to represent the two strands of
// genetic information in each genotype.
//

#ifndef _strand_h_
#define _strand_h_

#include <string>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "genome.h"


// A strand is a set of loci.  The implementation defined by
// this header file has the same interface as other implementations,
// but uses a simple vector to represent the loci in the strand.

class Locus;

class Strand {
friend class StrandIter;
public:
  Strand();
  ~Strand();
  static void set_chromosome_length(int n);
  fitness_t set_gene(int x, int y, mutation_t s);
  mutation_t get_gene(int x, int y);
  fitness_t get_fitness(int x);
  fitness_t strand_fitness();
  Strand* unravel(int n, int y);
  fitness_t merge(Strand *s1, Strand *s2);
  void print(std::ostream &sout, int limit);
private:
  Locus *v;
  int nl;				/* number of loci in strand */
  int ml;				/* size of v (max loci in this strand) */
  static int lmax;			/* max loci allowed in any strand */
  static int lseg;			/* allocation segment size */
  void extend_vector();			/* called if we need more loci */
  void move_down(int i);		/* make room for new locus */
  fitness_t copy_locus(int i, Locus *p, int j);
  fitness_t merge_locus(int i, Locus *p1, int j, Locus *p2, int k);
  Strand(const Strand &s);		/* disallow assignment and copying of strands */
  Strand operator=(const Strand &s);	/* these two operators are not implemented */
};

// An iterator for a Strand object returns a sequence of indices of
// mutated (non-zero) loci.  The indices are returned in order, from
// lowest to highest. Usage:
//	StrandIter si(Strand *sp)	create iterator si for strand sp
//	si				current index value [integer]
//	++si				increment to index of next non-zero locus
// The following code is an example of how to multiply together
// the values in a strand *sp:
//
//   for (StrandIter si(sp); si < max; ++si)
//     w *= sp.get_fitness(si);
//
// (where max is an externally defined value known to be higher than any
// locus index)

class StrandIter {
public:
  StrandIter(Strand &sp);
  int operator++();
  operator const int() const;
private:
  Locus *pv;
  int cur;
  int max;
};

#endif  // _strand_h_
