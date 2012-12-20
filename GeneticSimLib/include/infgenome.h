
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

// 
// infgenome.h -- "infinite genome" representation
//
// The implementation of add_mutations will guarantee that each
// new mutation goes into a locus that does not contain another
// mutation.  As the simulation advances, the genome effectively
// grows as long as it needs to be to accomodate new mutations,
// hence the name "infinite genome".
//
// In the actual implementation a "garbage collector", hidden from
// users, reclaims loci that are not being used.  The simulations
// will reach a steady state where new mutations simply replace
// older ones and the amount of space used is actually constant.
//

#ifndef _infgenome_h_
#define _infgenome_h_

#include "genome.h"

class ostream;			/* <iostream.h> */
class ISet;			/* "individual.h" */

class InfGenomeParamBlock : public GenomeParamBlock {
public:
  double d;			/* dominance factor */
};


typedef unsigned char *strand;

class InfiniteGenome : public Genome {
public:
  InfiniteGenome();
  ~InfiniteGenome();
  // virtual functions defined in Genome:
  fitness_t combine(Genome *p1, Genome *p2);
  fitness_t set_locus(int i, mutation_t s[]);
  fitness_t add_mutations(int n);
  fitness_t fitness();
  void print(ostream &s);
  Genome* operator =(Genome *x);
  // functions specific to this class (must be defined in all classes derived from Genome):
  static void reset_class();
  static void set_parameters(InfGenomeParamBlock *p);
  static void init_mutations(ISet &s);
  static void class_status(ostream &s);
  static void audit();
protected:
  // "virtual constructor"
  Genome *new_genome();
  // instance variables that define each infinited genome object
  strand s1;			/* pointer to strand 1 */
  strand s2;			/* pointer to strand 2 */
  fitness_t w;			/* relative fitness of these genes */
  InfiniteGenome *next;		/* 'next' and 'prev' implement a doubly linked */
  InfiniteGenome *prev;		/*    list of genome objects, used by garbage collection */
  // class variables
  static mutation_t s;		/* simulation parameter: selection coefficient */
  static double d;		/* simulation parameter: dominance */
  static double u;		/* simulation parameter: mutation rate */
  static int maxgl;		/* max genome length (bytes) */
  static int gl;		/* current genome length (bytes) */
  static strand free_loci;	/* free locus map */
  static int flp;		/* free locus pointer */
  static InfiniteGenome *head;	/* first genome in linked list */
  static int count;		/* number of genotypes currently active */
  static int garbage_collect();	/* garbage collection function */
  friend fitness_t fitness_function(int n1, int n2);
  friend fitness_t fitness_function(int n);
};

#endif

