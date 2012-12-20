
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

// Coho salmon simulation parameters.
//
// To add a new parameter, give it an ID and add it it to this definition of
// the enumeration type.  Add a line for the new parameter in the ParamTable
// structure so it can be read from the .rc file and/or the command line.  
//

typedef enum {
  NPopulations,			// number of populations to simulate
  MaxGenerations,		// max generations to create in any population
  PopSizeM,			// mean population size (for new generations)
  PopSizeSD,			// std dev of population size
  ReproRateM,			// mean reproductive rate
  ReproRateSD,			// std dev of reproductive rate
  FemaleP,			// percentage of newborns that are females
  MaleP,			// percentage of newborns that are males
  JackP,			// newborns that are jacks (males that return as 2-year olds)
  NSParams,			// last simulation parameter
  CeresNP,			// execution environment: number of jobs to
  PyrosNP,			//    run on ceres (Sun workstation) or pyros,
  AerosNP,			//    aeros, and geos (SGI trio).  Users can run
  GeosNP,			//    one job on ceres or, via MPI, up to 30 jobs
  NParams			//    at a time on the trio
} paramtype;

// Note: the following data structures are defined in both the front-end and
// back-end; since they are separate executables this file can be included
// by both....

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
    100.0, 0.0, MAXFLOAT },
  { MaxGenerations,
    "MaxGenerations",
    "-maxgen",
    "Max generations.................. ",
    100.0, 0.0, MAXFLOAT },
  { PopSizeM,       
    "PopSizeM",
    "-k",
    "Log population size.............. ", 
    3.47,   0.0, MAXFLOAT },
  { PopSizeSD,
    "PopSizeSD",
    "-sdk",
    "",
    0.0,   0.0, MAXFLOAT },
  { ReproRateM,
    "ReproRateM",
    "-r",
    "Log reproductive rate............ ",
    0.69,   0.0, MAXFLOAT },
  { ReproRateSD,
    "ReproRateSD",
    "-sdr",
    "",
    0.0,   0.0, MAXFLOAT },
  { FemaleP,
    "FemaleP",
    "-f",
    "Percent female................... ",
    0.401, 0.0, 1.0 },
  { MaleP,
    "MaleP",
    "-m",
    "Percent male..................... ",
    0.506, 0.0, 1.0 },
  { JackP,
    "JackP",
    "-j",
    "Percent jack..................... ",
    0.093, 0.0, 1.0 },
  { NSParams,"","","",0.0,0.0,0.0},
  { CeresNP,
    "CeresNP",
    "-ceres",
    "", 1, 0, 1 },
  { PyrosNP,
    "PyrosNP",
    "-pyros",
    "", 0, 0, 8 },
  { AerosNP,
    "AerosNP",
    "-aeros",
    "", 0, 0, 12 },
  { GeosNP,
    "GeosNP",
    "-geos",
    "", 0, 0, 10 }
};
