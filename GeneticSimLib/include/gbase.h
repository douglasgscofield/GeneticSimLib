
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

//
// gbase.h -- header file for genotype database operations
//
// A genotype database is a file with a large number (200,000 or more) 
// of genotypes which can be used to initialize individuals in a
// population.  The genes are stored in a compact representation in
// a file created by the 'mkgbase' application.  An object of type GBase
// reads a subset of the file and creates an internal database; the
// resulting object can be passed as a parameter to one of the
// Population class constructors to initialize the genes of individuals
// in the new population.  The intermediate step of creating a GBase
// object (instead of having the Population constructor read the file)
// makes it easier for parallel applications to read the file once and
// make copies of the GBase object for each process.
//
// A genotype database has the following components (both in the disk
// file and the GBase object):
//   A file header, containing the database ID and three integers that
//     give the number of items in the three segments of the database
//   The genome parameters used to create the database
//   An array s of mutation values, where s[i] is the effect of a
//     mutation on locus i.
//   A set of packed genotypes.  Each genotype in this section is 
//     a sequence of 16-bit words.  The first word is the number of 
//     mutated loci in the genotype; the remaining words are packed 
//     16-bit structs, one per locus, containing the index of the 
//     mutated locus and a flag that indicates which strand(s) 
//     contain the mutation.
//

#ifndef _gbase_h_
#define _gbase_h_

// Identifying tag expected as the first line of the genotype database:

#define GBTAG "CritterBase 1.0"

// The database has three segments.  The first three numbers in the
// file after the tag line are the sizes of each segment:

const int GB_NSEGS = 3;

typedef enum {
  GB_PARAM_SIZE,		/* number of genome parameters stored in file */
  GB_MUT_ARRAY_SIZE,		/* size of array of mutation effects */
  GB_SIZE			/* number of genotypes in the database */
} gb_segs;

// Genome parameters used to create the database.  This all-purpose
// set of definitions includes parameter names, default values, range
// of legal values, and command line names for the parameters, and can
// be used both by programs that make databases and programs that
// read values from it.

typedef float gb_param;

typedef enum {
  GB_NG,			// number of genotypes in a file
  GB_GL,			// genome length (number of loci)
  GB_U,				// gametic mutation rate
  GB_S,				// mean mutation effect
  GB_SDS,			// std dev of mutation effects
  GB_NPARAMS
} gb_param_type;

// The array of mutation effects is just an array of floats:

typedef float gb_mutation;

// A mutated locus is represented internally as a structure with a
// locus index and a strand indicator, and externally as a single
// 4-digit unsigned integer:

typedef enum {
  GB_HOMOZYGOUS_WILD,
  GB_HETEROZYGOUS_1,
  GB_HETEROZYGOUS_2,
  GB_HOMOZYGOUS_MUTANT
} gb_locus_type;

typedef union {
  struct {
    unsigned int i : 14;
    unsigned int f : 2;
  } s;
  unsigned short x;
} gb_locus;
  
// The internal representation of a genotype:

typedef struct {
  int n;
  gb_locus *m;
} gb_genotype;

// A GBase object has the same basic structure as the genotype
// database file.  Once it is created, the member functions can be
// used to access a specified individual and loci within the
// current individual.

class GBase {
public:
  GBase();
  GBase(int n, char *gfn);
  ~GBase();
  int size();			/* returns number of genotypes */
  void select(int n);		/* make genotype n the current one */
  int num_mutations();		/* number of mutated loci in current genotype */
  int next_locus();		/* bump index to next mutated locus */
  int locus_index();		/* index of current mutation */
  int locus_type();		/* heterozygous or homozygous */
  mutation_t locus_val(int n);	/* mutation value on strand n of current locus */
  mutation_t ms();		/* accessor for mean s used to build database */
  mutation_t sds();		/* ditto std deviation s */
  double u();			/* ditto mutation rate */
  int gl();			/* ditto genome length */
protected:
  gb_param *params;		/* genome parameters used to build database */
  gb_mutation *sv;		/* sv[i] is effect of mutation at locus i */
  gb_genotype *gv;		/* gv[i] is genotype i */
  int gb_size;			/* number of genotypes in the database */
  int cgx;			/* current genotype index */
  int clx;			/* current locus index */
};

#endif  // _gbase_h_
