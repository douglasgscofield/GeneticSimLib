
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

/*
  Coho salmon mutation accumulation simulation 

  Sequential version with command line interface.
*/

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>

#include "gsl.h"
#include "coho.h"
#include "cparams.h"		// simulation parameter definitions


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

void read_rcfile(double pb[]) {
  char line[100];
  std::ifstream rcfile(".cmasrc");
  double x;
  char *pn, *val;
  int i, valid;
  static char whitespace[] = " \t\n";
  static const char tagline[] = "# cmas 2.0";

  if (rcfile.bad())
    return;			// file doesn't exist? just use defaults

  rcfile.getline(line,sizeof(line));
  if (strcmp(line,tagline)) {
    std::cerr << "Corrupt .cmasrc file?  First line must be '" << tagline << "'" << std::endl;
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
      if (strcmp(ParamTable[i].name.c_str(),pn) == 0) {
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
	std::cerr << "Bad parameter value: " << pn << ", " << val << std::endl;
    else
      std::cerr << "Unknown parameter file entry: " << pn << "; ignored" << std::endl;
  }
}

void read_cmnd_line(int argc, char *argv[], double pb[]) {
  int i, j;
  double x;
  int valid;

  for (i = 1; i < argc; i++) {
    for (j = 0; j < NParams; j++)
      if (strcmp(argv[i],ParamTable[j].cmnd.c_str()) == 0) {
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
	    std::cerr << "Bad parameter value: " << argv[i-1] << ", " << argv[i] << std::endl;
    else
      std::cerr << "Unknown command: " << argv[i] << "; ignored" << std::endl;
  }
}

void load_parameters(int argc, char *argv[], double pb[]) {
  int i;

  for (i = 0; i < NParams; i++)
    pb[ParamTable[i].key] = ParamTable[i].val;
  read_rcfile(pb);
  read_cmnd_line(argc,argv,pb);
}

// print_parameters(double pb[]) -- print the parameters on the output
// stream so users can keep a record of experiment settings

void print_param(double pb[], paramtype p, int cr) {
  int i;

  for (i = 0; i < NParams; i++)
    if (ParamTable[i].key == p)
      break;

  if (cr)
    std::cout << std::endl << ParamTable[i].desc << pb[p];
  else
    std::cout << " \261 " << pb[p];
}

void print_parameters(double pb[]) {
  std::cout << "CMAS 2.0" << std::endl;

  print_param(pb,NPopulations,1);
  print_param(pb,MaxGenerations,1);
  print_param(pb,PopSizeM,1);
  print_param(pb,PopSizeSD,0);
  print_param(pb,ReproRateM,1);
  print_param(pb,ReproRateSD,0);
  print_param(pb,FemaleP,1);
  print_param(pb,MaleP,1);
  print_param(pb,JackP,1);

  std::cout << std::endl << std::endl;
}

int AGSize[NYEARS-1];		// sizes of ancestral generations

// This is gross.  Figure out a better way to specify the 
// fitness distribution....

#define WDistFile "../../lib/wdist.dat"

void init(double pb[]) {
  int i, n;
  CDFRNG *r;
  VirtualGenomeParamBlock vpb;

  r = new CDFRNG(WDistFile);
  if (r->bad())
    throw COHO_INIT_ERROR;

  vpb.R = r;
  VirtualGenome::set_parameters(&vpb);
  VirtualGenome::class_status(std::cout);

  LogNormalLogRNG rk(pb[PopSizeM],pb[PopSizeSD]);

  n = 0;
  for (i = 0; i < NYEARS-1; i++) {
    AGSize[i] = int(++rk + 0.5) * exp(0.5);
    n += AGSize[i];
  }
}

// Using the current parameters in the vector 'params' create a new population
// object for simulation

CohoPopulation *new_population(int pid, double *params) {
  CohoParamBlock pb;
  CohoPopulation *p;
  char pidstring[20];

  pb.kmax = params[PopSizeM];
  pb.sdk = params[PopSizeSD];
  pb.rmax = params[ReproRateM];
  pb.sdr = params[ReproRateSD];
  pb.u = 0.5;
  pb.pfemale = params[FemaleP];
  pb.pmale = params[MaleP];
  pb.pjack = params[JackP];

  sprintf(pidstring,"P%02d:%.0f",pid,params[NPopulations]);

  p = new CohoPopulation(&pb,AGSize,pidstring);

  return p;
}

void run(double *params) {
  CohoPopulation *p;
  ResultBlock *r;
  int maxgen = params[MaxGenerations];
  int i, n;
  double x;

  n = params[NPopulations];

  for (i=0; i<n; i++) {
    p = new_population(0,params);
    r = p->run(maxgen);
    delete p;
    delete r;
  }
}

int main(int argc, char *argv[]) {
  double params[NParams];
  FILE *pipe;

  load_parameters(argc,argv,params);
  print_parameters(params);

  init(params);

  try {
    run(params);
  }

  catch (COHO_ALLOC_ERROR e) {
    switch (e) {
    case COHO_INIT_ERROR:		std::cerr << "Initialization error" << std::endl; break;
    case COHO_NO_COHO:			std::cerr << "Couldn't allocate individual" << std::endl; break;
    case COHO_NO_NEW_GENERATION:	std::cerr << "Couldn't allocate generation" << std::endl; break;
    case COHO_NO_PARENTS:		std::cerr << "Couldn't allocate parent set" << std::endl; break;
    case COHO_NO_CHILD:			std::cerr << "Couldn't allocate new child" << std::endl; break;
    default:
      std::cerr << "Unknown runtime error" << std::endl;
    }
  }

  catch (...) {
    std::cerr << "Unknown (non-integer) runtime error" << std::endl;
  }
}
