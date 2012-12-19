/*
  Copyright (c) 1997 by the University of Oregon.
  See the COPYRIGHT file in this directory for permission
  to use and distribute this software.
*/

/*
  Mutation Accumulation Simulation

  Plain command line interface -- gets default settings from .rc file,
  then checks command line for further settings.  Prints results on
  stdout (cout).  
*/

#include <stdlib.h>
#include <ctype.h>
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>

#include "gsl.h"		// genetic simulation library

#define CDATE "7/97"

// Simulation parameters.  To add a new parameter, give it an ID and add it
// it to the enumeration type, then add a line for it in ParamTable, listing
// its ID, the string of its ID (the name used in the .rc file), the string
// used when specifying the value on the command line, the description
// used when it is printed for the user, and default, minimum, and maximum values.

typedef enum {
  NPopulations,			// number of populations to simulate
  MaxGenerations,		// max generations to create in any population
  PopSizeM,			// mean population size (for new generations)
  PopSizeSD,			// std dev of population size
  ReproRateM,			// mean reproductive rate
  ReproRateSD,			// std dev of reproductive rate
  MutEffectM,			// mean mutation effect
  MutEffectSD,			// std dev of mutation effect
  MutationRate,			// avg number of new mutations per individual
  GenomeLength,			// total number of loci in genome
  NChromosomes,			// number of chromosome pairs
  NCrossovers,			// avg number of cross-overs per combination
  Interval,			// used to print feedback on progress of runs
  NParams
} paramtype;

typedef struct {
  paramtype key;
  char *name;
  char *cmnd;
  char *desc;
  double val;
  double min;
  double max;
} param;

param ParamTable[] = {
  { NPopulations,   
    "NPopulations",
    "-pops",
    "Number of populations............ ", 
    250.0, 0.0, MAXFLOAT },
  { MaxGenerations,
    "MaxGenerations",
    "-maxgen",
    "Max generations.................. ",
    100.0, 0.0, MAXFLOAT },
  { PopSizeM,       
    "PopSizeM",
    "-k",
    "Population size.................. ", 
    2.0,   0.0, MAXFLOAT },
  { PopSizeSD,
    "PopSizeSD",
    "-sdk",
    "",
    0.0,   0.0, MAXFLOAT },
  { ReproRateM,
    "ReproRateM",
    "-r",
    "Reproductive rate................ ",
    2.0,   0.0, MAXFLOAT },
  { ReproRateSD,
    "ReproRateSD",
    "-sdr",
    "",
    0.0,   0.0, MAXFLOAT },
  { MutationRate,
    "MutationRate",
    "-u",
    "Mutation rate.................... ",
    0.5,   0.0, MAXFLOAT },
  { MutEffectM,
    "MutEffectM",
    "-s",
    "Mutation effect.................. ",
    0.05,  0.0, 1.0 },
  { MutEffectSD,
    "MutEffectSD",
    "-sds",
    "",
    0.02,  0.0, 1.0 },
  { GenomeLength,
    "GenomeLength",
    "-gl",
    "Genome Length.................... ",
    10000.0, 1.0, MAXFLOAT },
  { NChromosomes,
    "NChromosomes",
    "-nchrom",
    "Number of Chromosomes............ ",
    1.0, 1.0, MAXFLOAT },
  { NCrossovers,
    "NCrossovers",
    "-nx",
    "Mean number of cross-overs....... ",
    3.0, 0.0, MAXFLOAT },
  { Interval,
    "Interval",
    "-stars",
    "",
    10.0, 0.0, MAXFLOAT }
};

// load_parameters(double *pb) -- fill the array pb with the values of the
// simulation parameters.  Values are set in the following order:
//   (1) Initialize the block with default values.  
//   (2) If it exists, read the .rc file (procedure read_rcfile())
//   (3) Use any command line arguments (procedure read_cmnd_line())
// Note that this order gives command line arguments highest precedence.

int numberp(char *s) {
  int i;

  i = 0;
  while (isdigit(s[i]) || (s[i] == '.'))
    i += 1;

  return (i == strlen(s));
}

void read_rcfile(double *pb) {
  char line[100];
  ifstream rcfile(".masrc");
  double x;
  char *pn, *val;
  int i, valid;
  static char whitespace[] = " \t\n";

  if (rcfile.bad())
    return;			// file doesn't exist? just use defaults

  rcfile.getline(line,sizeof(line));
  if (strcmp(line,"# mas 1.0")) {
    cerr << "Corrupt .masrc file?  First line must be '#mas 1.0'" << endl;
    return;
  }

  while (rcfile.getline(line,sizeof(line))) {
    if (line[0] == 0)
      continue;
    if (line[0] == '#')
      continue;

    pn = strtok(line,whitespace);
    if (pn)
      val = strtok(NULL,whitespace);

    for (i = 0; i < NParams; i++)
      if (strcmp(ParamTable[i].name,pn) == 0) {
	valid = 1;
	if (!numberp(val))
	  valid = 0;
	else {
	  x = atof(val);
	  if ((x < ParamTable[i].min) || (x > ParamTable[i].max))
	    valid = 0;
	}
	break;
      }
    if (i < NParams)
      if (valid)
	pb[ParamTable[i].key] = x;
      else
	cerr << "Bad parameter value: " << pn << ", " << val << endl;
    else
      cerr << "Unknown parameter file entry: " << pn << "; ignored" << endl;
  }
}

void read_cmnd_line(double *pb, int argc, char *argv[]) {
  int i, j;
  double x;
  int valid;

  for (i = 1; i < argc; i++) {
    for (j = 0; j < NParams; j++)
      if (strcmp(argv[i],ParamTable[j].cmnd) == 0) {
	i = i+1;
	valid = 1;
	if (!numberp(argv[i]))
	  valid = 0;
	else {
	  x = atof(argv[i]);
	  if ((x < ParamTable[j].min) || (x > ParamTable[j].max))
	    valid = 0;
	}
	break;
      }
    if (j < NParams)
      if (valid)
	pb[ParamTable[j].key] = x;
      else
	cerr << "Bad parameter value: " << argv[i-1] << ", " << argv[i] << endl;
    else
      cerr << "Unknown command: " << argv[i] << "; ignored" << endl;
  }
}

void load_parameters(double *pb, int argc, char *argv[]) {
  int i;

  for (i = 0; i < NParams; i++)
    pb[ParamTable[i].key] = ParamTable[i].val;
  read_rcfile(pb);
  read_cmnd_line(pb,argc,argv);
}

// print_parameters(double pb) -- print the parameters on the output
// stream so users can keep a record of experiment settings

void print_param(double *pb, paramtype p, int cr) {
  int i;

  for (i = 0; i < NParams; i++)
    if (ParamTable[i].key == p)
      break;

  if (cr)
    cout << endl << ParamTable[i].desc << pb[p];
  else
    cout << " \261 " << pb[p];
}

void print_parameters(double *pb) {
  int i;

  cout << "MAS 3.0 [" << CDATE << "]" << endl;

  print_param(pb,NPopulations,1);
  print_param(pb,MaxGenerations,1);
  print_param(pb,PopSizeM,1);
  print_param(pb,PopSizeSD,0);
  print_param(pb,ReproRateM,1);
  print_param(pb,ReproRateSD,0);
  print_param(pb,MutEffectM,1);
  print_param(pb,MutEffectSD,0);
  print_param(pb,MutationRate,1);
  print_param(pb,GenomeLength,1);
  print_param(pb,NChromosomes,1);
  print_param(pb,NCrossovers,1);

  cout << endl << endl;
}

// init() -- miscellaneous initializations.
//  (1) initialize the genome class

void init(double *pb) {
  SparseGenomeParamBlock spb;

  spb.u = pb[MutationRate]; 
  spb.s = pb[MutEffectM];
  spb.sds = pb[MutEffectSD];
  spb.gl = 10000;
  spb.nchromosomes = 1;
  spb.maplength = 3;
  SparseGenome::set_parameters(&spb);
}

// Using the current parameters in the vector 'params' create a new population
// object for simulation

Population *new_population(double *params) {
  PopulationParamBlock pb;
  Population *p;

  pb.kmax = params[PopSizeM];
  pb.sdk = params[PopSizeSD];
  pb.rmax = params[ReproRateM];
  pb.sdr = params[ReproRateSD];
  pb.u = params[MutationRate];
  p = new Population(&pb);

  return p;
}

void run(double *params) {
  ResultBlock *r;
  Statistics s;
  int i, n, maxgen;
  Population *the_population;
  int stars = params[Interval];

  n = params[NPopulations];
  maxgen = params[MaxGenerations];

  for (i=0; i<n; i++) {
    the_population = new_population(params);
    r = the_population->run(maxgen);
    s.save(r->ngen);
    delete the_population;
    delete r;
    if (((i+1) % stars) == 0)
      cerr << "*";
    if (((i+1) % (stars * 80)) == 0)
      cerr << endl;
  }

  cerr << endl << endl;

  StatsBlock *sr = s.get_results();
  cout << "N................................ " << sr->n << endl;
  cout << "Mean extinction time (te)........ " << sr->mean << endl;
  cout << "Standard deviation of te......... " << sr->sd << endl;
  cout << "Coefficient of variation ........ " << sr->cv << endl;
  cout << "Min generations to extinction.... " << sr->min << endl;
  cout << "Max generations to extinction.... " << sr->max << endl;
  cout << endl;
}

int main(int argc, char *argv[]) {
  double params[NParams];

  load_parameters(params,argc,argv);
  print_parameters(params);
  init(params);
  run(params);
}
