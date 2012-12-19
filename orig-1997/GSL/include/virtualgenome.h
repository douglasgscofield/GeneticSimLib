
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

// 
// virtualgenome.h -- "virtual genome" representation
//
// This implementation does not store any genetic information at
// all.  An individual''s fitness is set by choosing a random
// number from an RNG object passed to the class initialization
// procedure.
//

#ifndef _virtualgenome_h_
#define _virtualgenome_h_

#include "genome.h"

class ostream;			/* <iostream.h> */
class ISet;			/* "individual.h" */
class RNG;			/* "rng.h" */

class VirtualGenomeParamBlock : public GenomeParamBlock {
public:
  RNG *R;
};

class VirtualGenome : public Genome {
public:
  VirtualGenome();
  ~VirtualGenome();
  // virtual functions defined in Genome:
  fitness_t combine(Genome *p1, Genome *p2);
  fitness_t set_locus(int i, mutation_t s[]);
  fitness_t add_mutations(int n);
  fitness_t fitness();
  void print(ostream &s);
  Genome* operator =(Genome *x);
  // functions specific to this class (must be defined in all classes derived from Genome):
  static void reset_class();
  static void set_parameters(VirtualGenomeParamBlock *p);
  static void init_mutations(ISet &s);
  static void class_status(ostream &s);
protected:
  Genome *new_genome();		// "virtual constructor"
  fitness_t w;			// fitness of this individual
  static RNG *R;		// RNG used to draw fitness values
};

#endif

