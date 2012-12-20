
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

//
// generation.cpp -- implementation of the Generation class.
//
// A generation is a collection of individuals.  Operations insert
// a new individual, return a pointer to an individual, remove an
// individual, and return the current size.
//
// In the basic implementation this class is little more than a wrapper
// for an ISet (set of Individuals).  One difference between an ISet and
// a Generation is in the destructor: when a Generation is deallocated 
// all the Individuals in that object are deallocated with it.  
//

#include <stddef.h>
#include <iostream>
#include "generation.h"
#include "individual.h"

// A new generation has no individuals; they'll be filled in by caller.

Generation::Generation() {
  is = new ISet;
}

// Deallocate the set and all the individuals in it:

Generation::~Generation() {
  Individual *pi;

  while ((pi = is->remove(0))) {
    delete pi;
  }

  delete is;
}

// Class initialization -- nothing to do ....

void Generation::reset_class() {
}

void Generation::set_parameters(GenerationParamBlock *pb) {
}

// Add x to the generation.

void Generation::insert(Individual *x) {
  is->insert(x);
}

// Remove the specified individual from the generation:

Individual *Generation::remove(int n) {
  return is->remove(n);
}

// Removes a random individual from the generation:


Individual *Generation::select() {
  return is->select();
}

// Index operator -- return a pointer to the nth individual:

Individual *Generation::operator [](int n) {
  return (*is)[n];
};

// simple accessor function....

int Generation::size() {
  return is->size();
}



