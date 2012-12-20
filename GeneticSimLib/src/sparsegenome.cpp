
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

// 
// sparsegenome.cpp -- implementation of the "sparse vector" representation
// of a genome.
//
// Linked list version: a locus contains an integer index and a float
// mutation value.  A strand is a linked list of loci, sorted according
// to index values.  A genome object consists of two strands, s1 and s2,
// and a float that stores the overall fitness (to avoid recomputation).
//
// Potential improvement: "smart copy" that represents a child as a set
// of new mutations plus pointers to parent genes; if the child won't
// be referenced until the parents are deallocated, use the parent
// destructor to reuse or copy loci as needed.
//

#include <string>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "strand.h"
#include "sparsegenome.h"
#include "individual.h"
#include "rng.h"

// Initialize class parameters with their default values

const mutation_t sdefault = 0.05;
const mutation_t sdsdefault = 0.0;
const double udefault = 1.0;	
const int gldefault = 1000;
const int nchromosomesdefault = 2;
const int chrlengthdefault = ceil(double(gldefault)/double(nchromosomesdefault));
const double maplengthdefault = 3.0;

mutation_t SparseGenome::s = sdefault;
mutation_t SparseGenome::sds = sdsdefault;
double SparseGenome::u = udefault;
int SparseGenome::gl = gldefault;
int SparseGenome::nchromosomes = nchromosomesdefault;
int SparseGenome::chrlength = chrlengthdefault;
double SparseGenome::maplength = maplengthdefault;

// This module uses two random number generators: one for mutation effects
// and one for the number of cross-overs in gene recombination.  The first
// is optional -- if it's not defined all mutations will have the same
// effect.

RNG *SparseGenome::rs = NULL;
RNG *SparseGenome::rx = NULL;

// Virtual constructor -- users call "make_genome()", defined in the base
// class.  The call is forwarded to this private member function to actually
// make a genome object.

Genome *SparseGenome::new_genome() {
  return (Genome *)new SparseGenome;
}

// The constructor creates a perfectly healthy individual (relative
// fitness = 1.0).  The constructor for the Strand class will create
// empty strands (no mutations).  

SparseGenome::SparseGenome() {
  sa = new Strand[nchromosomes];
  w = 1.0;
}

SparseGenome::~SparseGenome() {
  delete[] sa;
}

//  combine(x,y) -- combine the genes from genomes x and y and put
//  them in this genome.  Calculates the fitness of the individual
//  as it goes; stores and returns the fitness.

fitness_t SparseGenome::combine(Genome *p1, Genome *p2){
  SparseGenome *sp1 = (SparseGenome *)p1;
  SparseGenome *sp2 = (SparseGenome *)p2;
  fitness_t xw = 1.0;
  int i;
  int nx1, nx2;			// number of cross-overs from parent 1, parent 2
  Strand *lp1, *lp2;		// strands created by combining genes from p1, p2

  for (i = 0; i < nchromosomes; i++) {
    nx1 = ++(*rx);
    nx2 = ++(*rx);
    //    std::cout << nx1 << " " << nx2 << std::endl;
    //    sp1->sa[i].print(std::cout,10000);
    //    sp2->sa[i].print(std::cout,10000);
    lp1 = sp1->sa[i].unravel(nx1,0);
    lp2 = sp2->sa[i].unravel(nx2,1);
    xw *= sa[i].merge(lp1,lp2);
    //    sa[i].print(std::cout,10000);
    //    std::cout << sa[i].strand_fitness() << std::endl;
  }

  return w = xw;
}

// Assign the mutation values in s[0] and s[1] to the genotype at
// locus n.

fitness_t SparseGenome::set_locus(int n, mutation_t xs[]) {
  int i = n / chrlength;
  int j = n % chrlength;
  fitness_t xw = 1.0;

  for (int k = 0; k < 2; k++)
    if (xs[k])
      xw *= sa[i].set_gene(j,k,xs[k]);

  return w *= xw;
}

// add_mutations(int nm) -- add nm new mutations to this genome.
// get the values of these new mutations by calling new_mutation_value()
// Inline function new_mutation_value() generates a new mutation value; 
// if the effects are fixed, the value of the new mutation is s, the
// global selection coefficient.  Otherwise the effect is drawn from
// a random distribution.

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

fitness_t SparseGenome::add_mutations(int nm){
  int i, j, k, l, n;
  fitness_t xw = 1.0;

  for (i=0; i<nm; i++) {
    n = rword(gl);
    j = n / chrlength;
    k = n % chrlength;
    l = rword(2);
    xw *= sa[j].set_gene(k,l,new_mutation_value(s,rs));
  }
  return w *= xw;
}

// fitness() -- return the fitness value last calculated by add_mutations()
// or combine().

fitness_t SparseGenome::fitness() {
  return w;
}


// operator = is a deep copy operator; copy all the information
// in genome *x to this genome and return a pointer to this genome
// (which allows assignments to cascade).

Genome* SparseGenome::operator =(Genome *x){
  int i;
  SparseGenome *x1 = (SparseGenome *)x;

  //  *s1 = *(x1->s1);
  //  *s2 = *(x1->s2);
  //  w = x1->w;
  std::cerr << "sparse genome deep copy not implemented yet" << std::endl;
  return this;
}

  
// operator <<: put a written representation of a genome object into an 
// iostream (useful in demos and debugging).

void SparseGenome::print(std::ostream &sout) {
  static const int cutoff = 10;			// max loci to print
  int i;

  sout.setf(std::ios::fixed, std::ios::floatfield);
  sout.precision(4);

  sout << "sparse genome @ " << (void *)this << " w = " << w << std::endl << std::endl;
  for (i=0; i<nchromosomes; i++) {
    sout << "sa[" << i << "]: ";
    sa[i].print(sout,cutoff);
    sout << std::endl;
  }
}

// for debugging -- print out fitness of heterogeneous loci for a 
// few values:

void print_locus_fitnesses() {
  int i, j;
  mutation_t x, y;

  for (i=0; i<10; i++)
    for (j=0; j<=i; j++) {
      x = i*0.01;
      y = j*0.01;
      std::cout << "s( " << x << "," << y << ") = " << locus_fitness(x,y) << std::endl;
    }
}

// SparseGenome class static functions -- 
//   reset_class() reinitializes class variables
//   set_parameters() defines the operating parameters (e.g. the mutation
//     effects and mutation rate) and (re)initializes internal housekeeping
//     variables.  Its most important function is to create the exemplar object.
//   init_mutations() adds ancestral mutations to a set of individuals
//   class_status() is used in debugging and demo programs to print the
//     current state of the class

void SparseGenome::reset_class() {
  if (exemplar)
    delete exemplar;
  exemplar = (Genome *)new SparseGenome;
}

void SparseGenome::set_parameters(SparseGenomeParamBlock *p){
  if (p != NULL) {
    s = p->s;
    sds = p->sds;
    u = p->u * 2.0;			// convert gametic rate to zygotic rate
    gl = p->gl;
    nchromosomes = p->nchromosomes;
    chrlength = ceil(double(gl)/double(nchromosomes));
    maplength = p->maplength;
  }
  else {
    s = sdefault;
    sds = sdsdefault;
    u = udefault;
    gl = gldefault;
    nchromosomes = nchromosomesdefault;
    chrlength = chrlengthdefault;
    maplength = maplengthdefault;
  }

  Strand::set_chromosome_length(chrlength);

  if (rs)
    delete rs;
  if (sds > 0.0)
    rs = new GammaRNG(s,sds);
  else
    rs = NULL;
  
  if (rx)
    delete rx;
  rx = new PoissonRNG(maplength/nchromosomes);

  reset_class();
}

// init_mutations() initializes a set of individuals with "ancestral" mutations.

typedef Individual *IP;

void SparseGenome::init_mutations(ISet &is) {
  int p = 2*is.size();		// number of "strands" = 2 * #individuals
  int i, j, k, l;		// counters and temp indices
  int n;			// number of individuals in set
  mutation_t si;		// mutation effect of a single gene
  fitness_t xw;			// change in fitness due to this gene
  double qi;			// frequency of gene
  double gu;			// genic mutation rate
  int x, nx;			// id of individual, #individuals
  int nl = 0;			// for debugging, keep track of #loci and
  int nm = 0;			// #mutations affected by this step
  int *bin;			// also keep track of distribution of mutation values
  IP *v;
  SparseGenome *sgp;

  n = is.size();
  v = new IP[n];
  for (i=0; i<n; i++)		// copy individuals to local vector for faster access
    v[i] = is.remove(0);

  bin = new int[p+1];
  //  std::cout << "Placing initial mutations in " << p << " strands...." << std::endl;

  for (i=0; i<=p; i++)
    bin[i] = 0;

  gu = u/(2*gl);		// genic rate = individual rate / number of genes

  //  std::cout << "Genic mutation rate (U/2N) = " << gu << std::endl;

  // Iterate over all loci; pick a mutation effect for that locus, figure
  // out the expected frequency of a gene with that effect, and then add it
  // to selected individuals.  In this loop 'i' is the locus index,  'j'
  // is the chromosome number for that locus, and 'k' is the index of the locus
  // within the chromosome.

  // optimization: knowing that get_gene() and set_gene() in Strand scan
  // from "left to right", start with the higher locus indices and count
  // down to 0 -- this way the lookups and insertions will all be done
  // in one step....

  for (i = gl-1; i >= 0; i--) {
    j = i / chrlength;
    k = i % chrlength;
    si = new_mutation_value(s,rs); 	// get mutation effect
    qi = gu/hs(si);			// frequency of allele with effect si
    if (qi > 1.0)
      qi = 1.0;
    nx = int(rbinomial(qi,p));   	// number of genes expected to have this mutation
   
    if (nx > p)
      nx = p;

    //    if (qi) {
    //      std::cout << "locus " << i << ": ";
    //      std::cout << ", si = " << si;
    //      std::cout << ", hs = " << hs(si,0.0);
    //      std::cout << ", qi = " << qi;
    //      std::cout << ", nx = " << nx << std::endl;
    //    }

    nm += nx;				// record info on this mutation
    bin[nx] += 1;
    if (nx)
      nl += 1;

    while (nx) {
      x = rword(n);
      l = rword(2);
      sgp = (SparseGenome *)v[x]->genes;
      if (sgp->sa[j].get_gene(k,l) == 0.0) {
	xw = sgp->sa[j].set_gene(k,l,si);
	sgp->w *= xw;
	nx--;
      }
    }
  }

  // std::cout << "Distributed " << nm << " mutations across " << nl << " loci" << std::endl;
  // for (i=0; i<=p; i++)
  //   std::cout << i << ": " << bin[i] << std::endl;
}

// class_status prints the state of the class (no static members in this class,
// so nothing to print)

void SparseGenome::class_status(std::ostream &sout) {
  sout << "[Sparse Genome]" << std::endl;
  sout << "s = " << s << " +/- " << sds << "; u = " << u << "; ";
  sout << "G = " << gl << "; N = " << nchromosomes << "; ";
  sout << "L = " << chrlength << "; M = " << maplength << std::endl;    
}

