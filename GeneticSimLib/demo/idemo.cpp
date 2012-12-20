
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

// idemo -- demo for Individual and Generation classes; a Generation
// object uses the ISet (set of individuals) container class defined
// in individual.h.

#include <iostream>

#include "gsl.h"

const int gsize = 10;

void init();

int main(int argc, char *argv[]) {
  Generation moms;
  Generation dads;
  Generation x;
  Individual *mom;
  Individual *dad;
  Individual *child;
  UniformRNG r(0.0,1.0);
  int i;

  init();

  // Create 10 males and 10 females, and give them each 2 mutations.
  // Then create a child and put it in generation x:

  for (i=0; i<gsize; i++) {
    dad = new Individual;
    dad->set_sex(MALE);
    dad->genes->add_mutations(2);
    cout << dad->get_id() << ":  " << *dad;
    dads.insert(dad);

    mom = new Individual;
    mom->set_sex(FEMALE);
    mom->genes->add_mutations(2);
    moms.insert(mom);
    cout << mom->get_id() << ":  " << *mom;

    child = new Individual;
    child->set_sex((++r < 0.5) ? MALE : FEMALE);
    child->genes->combine(dad->genes,mom->genes);
    x.insert(child);
    cout << child->get_id() << ":  " << *child;

    cout << endl;
  }

  // Test the remove function and error detection in the Generation class.

  // Print the first three, remove the second, and print again:

  cout << "First three individuals in x: \n  ";
  for (i=0; i<3; i++)
    cout << x[i]->get_id() << "  ";
  cout << endl;
  child = x.remove(1);
  cout << "individual removed: " << child->get_id() << endl;
  cout << "First three individuals remaining in x: \n  ";
  for (i=0; i<3; i++)
    cout << x[i]->get_id() << "  ";
  cout << endl << endl;

  // Same test, but now remove the first one (special case):

  cout << "First three individuals in x: \n  ";
  for (i=0; i<3; i++)
    cout << x[i]->get_id() << "  ";
  cout << endl;
  child = x.remove(0);
  cout << "individual removed: " << child->get_id() << endl;
  cout << "First three individuals remaining in x: \n  ";
  for (i=0; i<3; i++)
    cout << x[i]->get_id() << "  ";
  cout << endl << endl;

  // There are now gsize-2 individuals in x.  Whittle it down to 2:

  for (i = 0; i < gsize-4; i++)
    child = x.remove(0);

  // Print the remaining two individuals:

  cout << "last two individuals in x: " << x[0]->get_id() << " and " << x[1]->get_id() << endl;

  // These should cause an error message:

  child = x[2];
  child = x.remove(2);

  // Test removing the last item:

  child = x.remove(1);
  cout << child->get_id() << endl;

  child = x.remove(0);
  cout << child->get_id() << endl;

  // Gen X should now be empty:

  cout << "Gen X size = " << x.size() << endl;
}

// Initializations -- set parameters used in genome class

void init() {
  SparseGenomeParamBlock pb;
  int gl = 10000;		// genome length (total #loci)
  int nchr = 1;			// number of diploid chromosomes

  pb.s = 0.05;
  pb.sds = 0.0;
  pb.u = 0.5;
  pb.gl = gl;
  pb.nchromosomes = nchr;
  pb.maplength = 3;
  SparseGenome::set_parameters(&pb);
}
