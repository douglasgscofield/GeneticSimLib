
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

// 
// genome.h -- representation of genes of a single individual
//
// class Genome is an abstract class that defines the basic
// operations on genes.  Applications should use these functions
// in order to be able to switch between representations.
//
// The class uses the "examplar" idiom ("Advanced C++ Programming
// Styles and Idioms" by J. C. Coplien) to implement a virtual
// constructor.  The use of the examplar gives users more
// flexibility in creating genome objects -- they can call the
// make_genome member function instead of naming the type of
// genome they want to create.
//

#ifndef _genome_h_
#define _genome_h_

#include <iostream>

class ISet;

// Use "mutation_t" in the definition of the effects of a single mutation.
// Use "fitness_t" for the cumulative effects of several mutations.  Both
// are single-precision floats to conserve space.

typedef double fitness_t;
typedef double mutation_t;

// The function hs(s) gives the relative reduction in fitness of
// a mutation when it is heterozygous with another mutation.  The
// parameter s is the relative fitness of the two genes.

inline fitness_t hs(mutation_t s) {
  return (0.5 * s) / (1 + 10.0 * s);
}

// This fitness function is used by several different genome
// implementations to compute the relative fitness of a single locus.

inline fitness_t locus_fitness(mutation_t s0, mutation_t s1) {
  mutation_t si = (s0 > s1) ? s0 : s1;
  mutation_t sj = (s0 > s1) ? s1 : s0;
  mutation_t s = (1 - (1-si)/(1-sj));

  return (1 - sj - hs(s));
}


// Every genome class uses at least the following parameters.  A derived
// class can defined its own parameter block as an extension of this class:

class GenomeParamBlock {
public:
  mutation_t s;			/* s = mean mutation effect */
  double u;			/* mu = gametic mutation rate */
};

class Genome {
protected:
  static Genome *exemplar;
  virtual Genome *new_genome() = 0;
public:
  friend Genome *make_genome();	/* "virtual constructor" */
  virtual ~Genome() {};
  virtual fitness_t combine(Genome *p1, Genome *p2) =0;
  virtual fitness_t set_locus(int i, mutation_t s[]) =0;
  virtual fitness_t add_mutations(int n) =0;
  virtual fitness_t fitness() =0;
  virtual void print(std::ostream &s) =0;
  virtual Genome* operator =(Genome *x) =0;
};

// In addition to the virtual functions listed above, which can be applied
// to any Genome object, the following functions must be defined in each
// derived class.  If there was such a thing as a virtual static function
// these would be defined that way -- functions that must be defined in each
// derived class, but applicable to the class as a whole and not any object.
//
//  static void set_parameters(GenomeParamBlock *p);
//     (sets params that are implementation-dependant)
//  static void init_mutations(ISet &s);
//     (puts "ancestral" mutations in objects in set s)
//  static void class_status(std::ostream &s);
//     (prints status of class on stream s)

// This operator is a "wrapper" that invoked the print member function:

inline std::ostream &operator <<(std::ostream &s, Genome *g) {
  g->print(s);
  return s;
}

#endif  // _genome_h_

