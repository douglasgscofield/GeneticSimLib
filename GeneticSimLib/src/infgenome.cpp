
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

// 
// infgenome.cpp -- implementation of the "infinite genome" representation
//
// Each locus is a pair of bits; an individual will have two arrays of
// bits.  The bits are stored as arrays of unsigned chars to make it easier
// to count mutations and do garbage collection.  
//


#include <string>
#include <stdlib.h>
#include <stddef.h>
#include <math.h>
#include <iostream>

#include "individual.h"		// needed for definition of ISet
#include "infgenome.h"
#include "rng.h"

// function prototypes of local functions

fitness_t fitness_function(int n1, int n2);
fitness_t fitness_function(int n);
strand build_free_locus_map(int maxgl);
unsigned short *build_bit_count();
char *strand_string(char *s, strand ps, int n);

// vars local to this implementation of the class
 
static unsigned short *bit_count =	// bit_count[i] = number of 1 bits in binary rep of i
   build_bit_count();

// Debugging -- audit number made, number deallocated

int num_allocated = 0;
int num_deallocated = 0;

// Virtual constructor -- users call "make_genome()", defined in the base
// class.  The call is forwarded to this private member function to actually
// make a genome object.

Genome *InfiniteGenome::new_genome() {
  return (Genome *)new InfiniteGenome;
}

// Constructor and destructor.  The garbage collector works by iterating
// over all current instances of a genome, so the constructor links new
// objects into a list of instances and the destructor removes the object.

InfiniteGenome::InfiniteGenome() {
  s1 = new unsigned char[maxgl];
  s2 = new unsigned char[maxgl];
  memset(s1,0,maxgl);
  memset(s2,0,maxgl);
  w = 1.0;
  prev = NULL;
  next = head;
  if (head)
    head->prev = this;
  head = this;
  count += 1;
  num_allocated += 1;
}
 
InfiniteGenome::~InfiniteGenome() {
  //  cout << "      deleting genome" << std::endl;
  delete s1;
  delete s2;
  if (prev)
    prev->next = next;
  if (next)
    next->prev = prev;
  if (head == this)
    head = next;
  count -= 1;
  num_deallocated += 1;
}

//  combine(x,y) -- combine the genes from genomes x and y and put
//  them in this genome.  Calculates the fitness of the individual
//  as it goes; stores and returns the fitness.

fitness_t InfiniteGenome::combine(Genome *p1, Genome *p2) {
  InfiniteGenome *ip1 = (InfiniteGenome *)p1;
  InfiniteGenome *ip2 = (InfiniteGenome *)p2;
  unsigned long rw;
  unsigned char r1;
  unsigned char r2;
  int i;
  int n1 = 0;
  int n2 = 0;
  unsigned char c1;
  unsigned char c2;

  for (i=0; i < gl; i++) {
    rw = rword();
    r1 = rw % 0x100;
    r2 = (rw >> 8) % 0x100;
    c1 = (((ip1->s1[i]) & r1) | ((ip1->s2[i]) & ~r1));
    c2 = (((ip2->s1[i]) & r2) | ((ip2->s2[i]) & ~r2));
    s1[i] = c1;
    s2[i] = c2;
    n1 += bit_count[c1 ^ c2];
    n2 += bit_count[c1 & c2];
  }

  w = fitness_function(n1,n2);
  return w;
}

fitness_t InfiniteGenome::set_locus(int i, mutation_t xs[]) {
    std::cerr << "InfiniteGenome::set_locus not implemented" << std::endl;
  return 0.0;
}

/*
  add_mutations(int nm) -- add nm new mutations to this genome.  

  The genome is organized as a vector of unsigned bytes.  The genome
  map, pointed to by 'free_loci', is also a vector of bytes.  If bit
  i in the map is 1, locus i is free to be used for a new mutation, 
  otherwise there is already a mutation at locus i in one or more 
  individuals.

  The class variable 'flp' is the index of the first byte that has
  a free locus.  To add new mutations, use the free locations in byte
  number flp.  If there are still mutations to add, increment flp and
  move on to the next byte.  This step might result in increasing
  the global gene length, 'gl'.

  The return value of the function is the new fitness of the individual.
  Since the new mutations are going to unused loci, they all have the
  effect of heterozygous mutations, so we can just update the current
  fitness (computed when the individual was "born") by (1-d*s)^nm.
*/

fitness_t InfiniteGenome::add_mutations(int nm) {
  unsigned char *flb;
  unsigned char *xb;
  unsigned char bits = 0;
  unsigned char mask;
  int n = 0;
  int i;

  w *= fitness_function(nm);			// true no matter where mutations inserted

  while (1) {
    flb = &free_loci[flp];
    xb = &s1[flp];

    if (nm <= bit_count[*flb]) { 		// if enough bits in current
      for (mask = 0x80; n < nm; mask >>= 1) 	// byte use them and return
	if (mask & *flb) {
	  bits |= mask;
	  n += 1;
	}
      *xb |= bits;
      *flb &= ~bits;

      return w;
    }

    *xb |= *flb;		// not enough; use all bits in current byte
    nm -= bit_count[*flb];	// and continue with the next byte
    *flb = 0;

    if (flp < gl-1) {
      flp += 1;			// still unused bits in current gene size
      continue;
    }
      
    // garbage collector will recycle unused bits, return index
    // of char with first unused bit (or -1 if no unused bits available)

    if ((i = garbage_collect()) >= 0) {
      flp = i;
      continue;
    }

    if (gl < maxgl) {
      flp += 1;			// extend gene and continue
      gl += 1;
      continue;
    }

    std::cerr << "genome::add_mutations: fatal error: out of loci\n";
    exit(1);
  }
}

// fitness() -- return the fitness value last calculated by add_mutations()
// or combine().

fitness_t InfiniteGenome::fitness() {
  return w;
}

// operator = is a deep copy operator; copy all the information
// in genome *x to this genome and return a pointer to this genome
// (which allows assignments to cascade).

Genome* InfiniteGenome::operator =(Genome *x) {
  //  return this;
    std::cerr << "deep copy not implemented in InfiniteGenome" << std::endl;
  return NULL;
}

// operator <<: put a written representation of a genome object into an iostream
// (useful in demos and debugging).  This version prints the genome on
// two lines; the application has to parse the string for the line
// separator (\n) if it wants to do anything fancy.
//
// This function and the class status function share a local constant named
// "indent" that tells how far to indent the bit vector.

static const std::string& indent = "        ";	// spaces to indent each strand
static const int linesize = 99;		// max cols to print in each strand
static const int nbytes = linesize/9; 	// print a byte and a space

void InfiniteGenome::print(std::ostream &sout) {
  char sbuf[linesize+2];
  const std::string& smindent = "  ";	// sigh -- set by hand to fill after fitness

  sout.setf(std::ios::fixed, std::ios::floatfield);
  sout.precision(4);

  sout << w << smindent << strand_string(sbuf,s1,(gl>nbytes)?(nbytes):(gl));
  sout << std::endl;
  sout << indent << strand_string(sbuf,s2,(gl>nbytes)?(nbytes):(gl));
  sout << std::endl;
}

// InfiniteGenome class static functions -- 
//   reset_class() reinitializes class variables
//   set_parameters() defines the operating parameters (e.g. the mutation
//     effects and mutation rate) and (re)initializes internal housekeeping
//     variables (e.g. garbage collection).  Its most important function,
//     though, is to create the exemplar object.
//   init_mutations() adds ancestral mutations to a set of individuals
//   class_status() is used in debugging and demo programs to print the
//     current state of the class

void InfiniteGenome::reset_class() {
  InfiniteGenome *pg;

  flp = 0;
  gl = 1;
  count = 0;

  if (free_loci)
    delete free_loci;
  free_loci = build_free_locus_map(maxgl);

  while (head) {
    pg = head;
    head = head->next;
    delete pg;
  }

  if (exemplar)
    delete exemplar;
  exemplar = (Genome *)new InfiniteGenome;
}

void InfiniteGenome::set_parameters(InfGenomeParamBlock *p) {
  if (p != NULL) {
    s = p->s / 2.0;		// new meaning of s: per locus effect
    d = p->d;			// no param input for dominance
    u = 2.0 * p->u;		// convert gametic to zygotic rate
  }
  else {			// no param block; use defaults
    s = 0.025;
    d = 1.0;
    u = 1.0;
  }

  reset_class();
}

// init_mutations(s) is used to initialize the set of individuals s with
// "ancestral" mutations. 

void InfiniteGenome::init_mutations(ISet &is) {
    std::cerr << "Initial mutations not defined for InfiniteGenome class" << std::endl;
}

void InfiniteGenome::class_status(std::ostream &sout) {
  int i, j;
  char sbuf[linesize+2];

  sout << "Inf Genome: ";
  sout << "count = " << count << ", ";
  sout << "gl = " << gl << ", ";
  sout << "flp = " << flp << ", ";
  sout << "s = " << s;
  sout << std::endl;

  sout << indent << strand_string(sbuf,free_loci,(gl>nbytes)?(nbytes):(gl));
  sout << std::endl;
}


//
// ------- Local Procedures ------------
//
// The following procedures are local to this implementation.  They
// help manage the class static variables, or are "helper" procedures
// for methods, etc.

//  Garbage collector -- return 1 if we can find some unused loci.

int InfiniteGenome::garbage_collect() {
  InfiniteGenome *p = head;
  int i;

  for (i = 0; i < gl; i++)
    free_loci[i] = 0xFF;

  while (p) {
    for (i = 0; i < gl; i++)
      free_loci[i] &= ~(p->s1[i] | p->s2[i]);
    p = p->next;
  }

  for (i = 0; i < gl; i++)
    if (free_loci[i])
      return i;

  return -1;
}


//
// This procedure initializes the static global bit_count so that
// bit_count[i] is the number of 1 bits in the binary representation
// of the integer i (0 <= i <= 255).
//

unsigned short audit_count[256];

unsigned short *build_bit_count() {
  int i, j, k;
  unsigned short *nb = new unsigned short[256];

  for (i = 0; i < 256; i++) {
    nb[i] = 0;
    for (j = 0, k = i; j < 8; j++, k >>= 1)
      if (k & 1)
	nb[i] += 1;
    audit_count[i] = nb[i];
  }

  return nb;
}

//
// Allocate the free locus map.  We want a vector of size genome_t[nwords],
// so it is as long as a genome, but it will be managed as a vector of chars
// so we can use bit_count and other utilities that are based on chars.
//

unsigned char *build_free_locus_map(int maxgl) {
  unsigned char *p =  new unsigned char[maxgl];
  int i;

  for (i=0; i<maxgl; i++)
    p[i] = 0xFF;

  return p;
}

// Build a printable string with a representation of a single strand
// of a genome.  p points to the string buffer to fill, ps points to
// the strand to print, and n is the number of segments (unsigned chars)
// to convert.  Return a pointer to the filled buffer (this lets the 
// string be sent directly to an output stream).

char *strand_string(char *p, strand ps, int n) {
  int i, j;
  unsigned char seg;
  char *ss = p;

  for (i=0; i<n; i++) {
    seg = ps[i];
    for (j=0; j<8; j++, seg <<= 1)
      if (seg & 0x80)
	*p++ = '1';
      else
	*p++ = '0';
    *p++ = ' ';
  }
  *p = '\0';

  return ss;
}
 
// Class (static) variables.  The first three can be given new values
// by calling init_infinite_genome_class.

mutation_t InfiniteGenome::s = 0.025;	// selection coefficient (effect of one mutation)
double InfiniteGenome::d = 1.0;		// dominance (relative effect of mutation over wild)
double InfiniteGenome::u = 1.0;		// zygotic mutation rate
int InfiniteGenome::maxgl = 1024;	// max genome length (bytes)

strand InfiniteGenome::free_loci = 	// 1 bits indicate loci where mutations can be added
   build_free_locus_map(InfiniteGenome::maxgl);
int InfiniteGenome::flp = 0;		// free locus pointer (index into free_loci)
int InfiniteGenome::gl = 1;		// current length (in bytes) of genome
InfiniteGenome *InfiniteGenome::head = NULL;		// linked list of allocated genomes
int InfiniteGenome::count = 0;			// number of genomes allocated

// Fitness function.  This version weights each mutation the same,
// with an effect defined by the simulation parameter 's'.  The fitness
// function is additive: heterozygous mutations have an effect of
// s, and homozygous mutations have an effect of 2s.  
//
// The fitness is thus
//	w = (1-d*s)^n1 * (1-2*s)^n2
// where n1 is the number of heterozygous loci (one mutation) and n2
// is the number of homozygous loci (two mutations).
//
// The parameter 'd' defines the effect of dominance.  If d = 1, the
// mutation is neutral, i.e. neither dominant nor recessive, and the
// contribution of this locus is 1-s.  If d = 0, the mutation is recessive,
// the the effect of the mutation is 0.  If d = 2, the mutation is
// dominant and should have full effect, as if the locus were homozygous,
// and thus the contribution of this locus is (1-2*s). 

fitness_t fitness_function(int n1, int n2) { 
  double d = InfiniteGenome::d;
  mutation_t s = InfiniteGenome::s;
  return (pow((1.0-d*s),n1) * pow((1.0-2.0*s),n2));
}

fitness_t fitness_function(int n) {
  double d = InfiniteGenome::d;
  mutation_t s = InfiniteGenome::s;
  return pow(1.0-d*s,n);
}


// Debugging aid ---

void InfiniteGenome::audit() {
  int i;

  std::cerr << "Audit: " << std::endl;

  std::cerr << "Checking global arrays..." << std::endl;
  for (i = 0; i < 256; i++) {
    if (audit_count[i] != bit_count[i]) {
        std::cerr << "error: bit_count[" << i << "] = " << bit_count[i] << std::endl;
    }
  }
  std::cerr << "done" << std::endl;

  std::cerr << "Counters: " << std::endl;
  std::cerr << "  objects in use:       " << count << std::endl;
  std::cerr << "  objects allocated:    " << num_allocated << std::endl;
  std::cerr << "  objects deallocated:  " << num_deallocated << std::endl;
}

