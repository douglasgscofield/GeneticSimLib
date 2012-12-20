
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

// 
// virtualgenome.C -- implementation of the "virtual genome" representation
//

#include <string.h>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <iostream.h>

#include "virtualgenome.h"
#include "rng.h"

// Class static variable: pointer to the random number generator
// that will draw fitness values

RNG *VirtualGenome::R = NULL;

// Virtual constructor -- users call "make_genome()", defined in the base
// class.  The call is forwarded to this private member function to actually
// make a genome object.

Genome *VirtualGenome::new_genome() {
  return (Genome *)new VirtualGenome;
}

// Constructor and destructor: nothing to do...

VirtualGenome::VirtualGenome() {
  w = 1.0;
}
 
VirtualGenome::~VirtualGenome() {
}

// combine(x,y) -- in this version, the fitness of the parents
// has no bearing on the fitness of the child -- just draw a
// random variable and use it as the fitness of this individual....

fitness_t VirtualGenome::combine(Genome *, Genome *) {
  w = ++(*R);
  return w;
}

// set_locus() should not be called by an application that uses
// a virtual genome:

fitness_t VirtualGenome::set_locus(int i, mutation_t xs[]) {
  cerr << "VirtualGenome::set_locus not implemented" << endl;
  return 0.0;
}

// Don't actually add mutations, but decrease the fitness by an
// equivalent amount

mutation_t new_mutation_value(mutation_t s, RNG *rs) {
  mutation_t x;
  if (rs) {
    x = ++(*rs);
    if (x > 1.0)
      return 1.0;
    else
      return x;
  }
  else
    return s;
}

fitness_t VirtualGenome::add_mutations(int nm) {
  fitness_t xw = 1.0;
  int i;
  static RNG *rs = new GammaRNG(0.05,0.02); // params used to create distribution...

  for (i = 0; i < nm; i++) 
    xw *= locus_fitness(new_mutation_value(0.05,rs),0.0);

  return w *= xw;
}

// fitness() -- return the fitness value last calculated by add_mutations()
// or combine().

fitness_t VirtualGenome::fitness() {
  return w;
}

// operator = is a deep copy operator; copy all the information
// in genome *x to this genome and return a pointer to this genome
// (which allows assignments to cascade).

Genome* VirtualGenome::operator =(Genome *x) {
  VirtualGenome *px = (VirtualGenome *)x;

  w = px->w;
  return this;
}

// The base class defines operator << which calls the print()
// member function defined separately in each derived class.
// For this class just print the fitness.  To be consistent with
// other implementations print an entire line.

void VirtualGenome::print(ostream &sout) {
  sout.setf(ios::fixed, ios::floatfield);
  sout.precision(4);
  sout << w << endl;
}

// VirtualGenome class static functions -- 

void VirtualGenome::reset_class() {
  if (exemplar)
    delete exemplar;
  exemplar = (Genome *)new VirtualGenome;
}

void VirtualGenome::set_parameters(VirtualGenomeParamBlock *p) {
  if (R)
    delete R;

  if (p) 
    R = p->R;
  else
    R = new NormalRNG(0.367,0.04);

  reset_class();
}

class ISet;

void VirtualGenome::init_mutations(ISet &is) {
  cerr << "Initial mutations not defined for VirtualGenome class" << endl;
}

void VirtualGenome::class_status(ostream &sout) {
  sout << "Virtual Genome: R = " << (void *)R << endl;
}
