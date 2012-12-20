
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

/*
  gdriver --- test driver for genome modules.

  In each iteration of the test:
  * create 'ng' new individuals by crossing genes from current individuals
  * add mutations to each new individual
  * deallocate the old generation and any new individuals with health
    less than the threshold
*/

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <strstream.h>
#include <iomanip.h>
#include <time.h>
#include <unistd.h>

#include "gsl.h"


// Define pointer to genome object:

typedef Genome *GenomeP;

// global vars

typedef enum {Infinite, Sparse} genome_type;
genome_type cur_genome_type = Infinite;

const int ng = 5;		// number of genome objects
const int reps = 10000;		// number of test runs to make

GenomeP g0[ng];			// arrays to hold generations
GenomeP g1[ng];
GenomeP *g;
GenomeP *og;

fitness_t threshold = 0.5;

clock_t start_time, end_time;
char *start_mem, *end_mem;

// local functions

void init(int argc, char *argv[]);
void new_generation_step();
void add_mutations_step();
void zap_individuals_step();

int main(int argc, char *argv[]) {
  int i;
  double t;

  init(argc,argv);		// also builds initial global vars

  start_time = clock();
  start_mem = (char *)sbrk(0);
  
  cout << "Starting..." << endl;

  for (i=0; i<reps; i++) {
    new_generation_step();
    add_mutations_step();
    zap_individuals_step();
  }
  cout << endl;

  end_time = clock();
  end_mem = (char *)sbrk(0);

  t = ((double)end_time - (double)start_time)/CLOCKS_PER_SEC;

  cout << "Done... " << endl;
  cout << "time = " << t << " seconds" << endl;
  cout << "memory: " << end_mem-start_mem << " bytes" << endl;

  return 0;
}

// Create the new generation -- do 'ng' random matings to fill the
// new generation with offspring from the current generation.

void new_generation_step() {
  int i;
  int x, y;
  fitness_t w, w0;

  GenomeP *t = g;		// swap the generation pointers
  g = og;
  og = t;

  w = 0.0;
  for (i=0; i<ng; i++) {
    g[i] = make_genome();
    x = rword(ng);
    y = rword(ng);
    w0 = g[i]->combine(og[x],og[y]);
    w += w0;
  }
}

// Add random mutations to each genome -- tests the add_mutations()
// member function of the genome class.

void add_mutations_step() {
  int i;
  int nm, n1, n2;
  static PoissonRNG rnm(1.0);

  for (i=0; i<ng; i++) {
    nm = ++rnm;
    g[i]->add_mutations(nm);
  }
}

// Delete all the individuals in the old generation, and replace
// any genome in the new generation that has a fitness below a set
// threshold.

void zap_individuals_step() {
  int i;
  int n1, n2;

  for (i=0; i<ng; i++)
    delete og[i];

  for (i=0; i<ng; i++)
    if (g[i]->fitness() < threshold) {
      delete g[i];
      g[i] = make_genome();
    }
}

void init(int argc, char *argv[]) {
  int i;

  if (argc > 1) {
    if (strcmp(argv[1],"-infinite") == 0)
      cur_genome_type = Infinite;
    else if (strcmp(argv[1],"-sparse") == 0)
      cur_genome_type = Sparse;
  }

  if (cur_genome_type == Infinite) {
    InfGenomeParamBlock ipb;
    ipb.s = 0.05;
    ipb.u = 0.5;
    ipb.d = 1.0;
    InfiniteGenome::set_parameters(&ipb);
  }
  else {
    SparseGenomeParamBlock spb;
    spb.u = 0.5;
    spb.s = 0.05;
    spb.sds = 0.02;
    spb.gl = 100;
    spb.nchromosomes = 1;
    spb.maplength = 3;
    SparseGenome::set_parameters(&spb);
  }

  g = g0;
  og = g1;

  for (i=0; i<ng; i++)
    g[i] = make_genome();
}


