
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

/*                                                                                  
  sparsegenome.h -- "sparse vector" representation of a genome                      
                                                                                    
  Terminology used in this module:                                                  
     genome:	   collection of all the genetic information in a species           
     genotype:     an instance of a genome, i.e. the genetic information in         
 		   a single individual from the population                          
     gene:	   abstract representation of functional portion of the genome;     
 		   represented by a real number between 0 and 1, which is the       
 		   relative fitness of the gene (1 = healthy).                      
     chromosome:   a sequence of genes; note that the genome may be broken          
 		   down into several chunks, so there may be several chromosomes    
 		   per individual                                                   
     haploid:	   one chromosome (or set of chromosomes) per individual            
     diploid:	   two chromosomes per individual, one inherited from each parent   
     locus:	   a specified location in the genome, consisting of the two genes  
 		   at that location in diploid populations; the relative fitness    
 		   of a locus is a function of the fitness of the two genes.        
 		                                                                    
  Implementation:                                                                   
     Genome:	   a set of diploid chromosomes, all of the same length             
     Strand:	   the data structure used to represent a chromosome, in this       
 		   case a linked list of loci                                       
     Locus:        a single object representing the two genes at a locus            
 		                                                                    
  Example:  If the genome has 1000 loci and 10 diploid chromosomes, it will be      
  represented in this program by 10 strands.  Each strand will be a linked list     
  of up to 100 loci each.  In the sparse vector representation, only loci that      
  have one or two mutated genes will be represented.                                
 		                                                                    
*/

#ifndef _sparsegenome_h_
#define _sparsegenome_h_

#include "genome.h"

class ostream;			/* <iostream.h> */
class ISet;			/* "individual.h" */
class RNG;			/* "rng.h" */

class Strand;

class SparseGenomeParamBlock : public GenomeParamBlock {
public:
  mutation_t sds;		/* std deviation of s (s = mutation effect) */
  int gl;			/* G = number of loci in genome */
  int nchromosomes;		/* N = number of haploid chromosomes */
  double maplength;		/* M = genetic map length (unit = Morgans) */
};

class SparseGenome : public Genome {
public:
  SparseGenome();
  ~SparseGenome();
  // virtual functions defined in Genome:
  fitness_t combine(Genome *p1, Genome *p2);
  fitness_t set_locus(int i, mutation_t s[]);
  fitness_t add_mutations(int n);
  fitness_t fitness();
  void print(ostream &s);
  Genome* operator =(Genome *x);
  // functions specific to this class (must be defined in all classes derived from Genome):
  static void reset_class();
  static void set_parameters(SparseGenomeParamBlock *p);
  static void class_status(ostream &s);
  static void init_mutations(ISet &is);
protected:
  // "virtual constructor"
  Genome *new_genome();
  // instance variables:
  Strand *sa;			/* array of strands */
  fitness_t w;			/* relative fitness of these genes */
  // class variables
  static mutation_t s;		/* mean mutation effect */
  static mutation_t sds;		/* std dev mutation effect */
  static double u;		/* mutation rate */
  static int gl;		/* total number of loci in genome */
  static int nchromosomes;	/* number of haploid chromosomes */
  static int chrlength;		/* number of loci per chromosome */
  static double maplength;	/* mean number of crossovers per genome */
  static RNG *rs;		/* RNG for variable mutation effects */
  static RNG *rx;		/* RNG for number of cross-overs per chromosome */
  friend mutation_t new_mutation_value();
};

#endif

