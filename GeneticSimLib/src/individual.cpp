
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

//
// Individual.cpp -- In this very basic implementation, an individual has
// a set of genes and a sex.  It also has a unique ID (or "tag") that is
// used mainly in debugging and demo programs.  Any more complex attributes 
// are left for derived classes to define.
//
// This module also defines an ISet, a container class for sets of individuals.
// Operations on ISets include inserting and removing individuals and obtaining
// the current size of the set.
//

#include <iostream>
#include <stddef.h>
#include "genome.h"
#include "infgenome.h"
#include "sparsegenome.h"
#include "individual.h"
#include "rng.h"
#include "gbase.h"

Genome *make_genome();

// This class variable is used to give each individual a unique ID:

int Individual::id = 0;

// Class initialization -- reset_class() and set_parameters() are 
// defined for consistency with other classes; there are no class
// variables to set...

void Individual::reset_class() {
}

void Individual::set_parameters(IndividualParamBlock *pb) {
}


Individual::Individual() : genes(make_genome()) {
  d_sex = NONE;
  d_id = id++;
}

Individual::~Individual() {
  delete genes;
}

void Individual::print(std::ostream &sout) {
  if (d_sex == NONE)
    sout << "[?]  ";
  else if (d_sex == MALE)
    sout << "[M]  ";
  else 
    sout << "[F]  ";
  sout << std::endl << genes;
}


void Individual::set_sex(sex_t s) {
  d_sex = s;
}

sex_t Individual::get_sex() {
  return d_sex;
}

int Individual::get_id() {
  return d_id;
}

// Copy the mutations in genotype 'n' in database 'gb' to this
// individuals genes:

void Individual::copy_genes(GBase *gb, int n) {
  int i, j, nm;
  mutation_t s[2];

  gb->select(n);
  nm = gb->num_mutations();

  for (i = 0; i < nm; i++) {
    j = gb->locus_index();
    s[0] = gb->locus_val(0);
    s[1] = gb->locus_val(1);
    genes->set_locus(j,s);
    gb->next_locus();
  }
}


// ------- ISet implementation ----------

// For debugging and tracking down memory leaks, keep track of the
// number of ISets allocated with this variable:

int ISet::count = 0;

// A set is a list tied together by "ILinks".  The ith individual
// is the ith element in the list.

ISet::ISet() {
  d_size = 0;
  d_first = NULL;
  count += 1;
}

// When a set is destroyed, only the "glue" that defines the set is
// deleted; the individual members of the set are NOT deleted.

ISet::~ISet() {
  ILink *pi = d_first;
  ILink *pt;

  while (pi) {
    pt = pi;
    pi = pi->next;
    delete pt;
  }

  count -= 1;
}

// Add x to the set.  Implements a simple "insert-front" insertion
// since there is no significance to an individual's place in the set.

void ISet::insert(Individual *x) {
  ILink *pi = new ILink;

  pi->p = x;			// allocate a link struct, fill it with
  pi->next = d_first;		// this object, and link it in to main list
  d_first = pi;

  d_size += 1;
}

// Remove the specified individual from the set:

Individual *ISet::remove(int n) {
  ILink *pi;
  ILink *pt;
  Individual *px;
  int i;

  if (n >= d_size)
    return NULL;

  if (n == 0) {
    pi = d_first;
    d_first = d_first->next;
  }
  else {
    i = 0;
    pt = d_first;
    while (i < n-1) {
      pt = pt->next;
      i += 1;
    }
    pi = pt->next;		// pt points to link before g[i]
    pt->next = pt->next->next;	// unsplice g[i]
  }

  // pi now points at the ILink for the individual to return; extract
  // the individual, destroy the link, and return

  px = pi->p;
  delete pi;

  d_size -= 1;
  return px;
}

// Select removes a random individual from the set:


Individual *ISet::select() {
  return remove(rword(d_size));
}

// Index operator -- return a pointer to the nth individual:

Individual *ISet::operator [](int n) {
  int i = 0;
  ILink *pi = d_first;

  if (n >= d_size)
    return NULL;

  while (i < n) {		// slow version for small sets: just
    pi = pi->next;		// iterate down the list....
    i += 1;
  }

  return pi->p;
};

// simple accessor function....

int ISet::size() {
  return d_size;
}
