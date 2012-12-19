
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

/*
  gdemo -- genome class demo

  This demo displays a small "population" and shows how new mutations
  are added.  It also demonstrates how the "garbage collection" method
  works in the case of variable length bit-vector genomes.  

  The program uses the XFORMS user interface library.
*/

#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <strstream.h>
#include <iomanip.h>

#include "gsl.h"

#include "forms.h"
#include "ggui.h"

FD_ggui *the_gui;

const int ng = 2;

Genome *g[ng];

// parameters and variables used in the genome display

const int linewidth = 100;
const int indent = 4;
const int bufsize = 500*linewidth;
const int statsize = 200;

char display_buffer[bufsize];

typedef enum {Infinite, Sparse} genome_type;
genome_type cur_genome_type = Infinite;

// local functions

void init_gui();
void init_obs();
void add_mutations_step();
void zap_individuals_step();
void genome_class_status(ostream &s);
void init_genome_class();

// Initialize the local objects (array of genomes) and user interface
// objects, then transfer control to the interface.  There are only
// two buttons:  'quit' and 'step'.  When the 'step' button is pushed,
// we will alternate between calling add_mutations_step(), which will
// add random new mutations to the genotypes, and zap_individuals_step(),
// which will kill off and replace some genotypes.

int main(int argc, char *argv[]) {
  fl_initialize(&argc,argv,"Genome Demo",0,0);

  if (argc > 1) {
    if (strcmp(argv[1],"-infinite") == 0)
      cur_genome_type = Infinite;
    else if (strcmp(argv[1],"-sparse") == 0)
      cur_genome_type = Sparse;
  }

  init_genome_class();
  init_obs();
  init_gui();

  while (1)
    fl_do_forms();

  return 0;
}

// Add random mutations to each genome -- tests the add_mutations()
// member function of the genome class.

void add_mutations_step() {
  int i;
  int nm, n1, n2;
  char num1[linewidth], num2[linewidth];
  ostrstream sout(display_buffer,bufsize);
  PoissonRNG rnm(1.0);

  fl_freeze_form(the_gui->ggui);
  fl_clear_browser(the_gui->Browser);

  genome_class_status(sout);

  for (i=0; i<ng; i++) {
    nm = ++rnm;
    sprintf(num1,"%3d ",nm);
    g[i]->add_mutations(nm);
    sout << g[i];
  }
  sout << ends;

  fl_addto_browser(the_gui->Browser,display_buffer);
  fl_unfreeze_form(the_gui->ggui);
}

// Destroy any genome that has a fitness below a set threshold and
// replace it with a new genome -- tests the fitness() member
// function of the genome class.

const fitness_t threshold = 0.5;

void zap_individuals_step() {
  int i;
  int n1, n2;
  char num[linewidth];
  ostrstream sout(display_buffer,bufsize);

  fl_freeze_form(the_gui->ggui);
  fl_clear_browser(the_gui->Browser);

  genome_class_status(sout);

  for (i=0; i<ng; i++) {
    if (g[i]->fitness() < threshold) {
      delete g[i];
      g[i] = make_genome();
    }
    sout << g[i];
  }
  sout << ends;

  fl_addto_browser(the_gui->Browser,display_buffer);
  fl_unfreeze_form(the_gui->ggui);
}

void genome_class_status(ostream &sout) {
  if (cur_genome_type == Infinite)
    InfiniteGenome::class_status(sout);
  else
    SparseGenome::class_status(sout);
  sout << endl;
}

void init_genome_class() {
  InfGenomeParamBlock ipb;
  SparseGenomeParamBlock spb;
  IndividualParamBlock pb;

  if (cur_genome_type == Infinite) {
    ipb.s = 0.05;
    ipb.u = 0.5;
    ipb.d = 1.0;
    InfiniteGenome::set_parameters(&ipb);
  }
  else {
    spb.s = 0.05;
    spb.u = 0.5;
    spb.sds = 0.02;
    spb.gl = 100;
    spb.nchromosomes = 4;
    spb.maplength = 3;
    SparseGenome::set_parameters(&spb);
  }

  Individual::set_parameters(&pb); 	// no parameters to fill....
}

void init_obs() {
  int i;
  Individual *ip;
  ISet is;

  for (i=0; i<ng; i++) {
    ip = new Individual;
    g[i] = ip->genes;
    is.insert(ip);
  }
}

void init_gui() {
  int i;
  ostrstream sout(display_buffer,bufsize);

  the_gui = create_form_ggui();

  // for some reason these two aren't being set by fdesign:
  fl_set_browser_fontsize(the_gui->Browser,FL_NORMAL_SIZE);
  fl_set_browser_fontstyle(the_gui->Browser,FL_FIXED_STYLE);  

  genome_class_status(sout);

  for (i=0; i<ng; i++)
    sout << g[i];
  sout << ends;

  fl_addto_browser(the_gui->Browser,display_buffer);
  fl_show_form(the_gui->ggui,FL_PLACE_CENTER,FL_FULLBORDER,"Genome Demo");
}

// adios....

void QuitButtonCB(FL_OBJECT *p0, long p1) {
  exit(0);
}

// Alternate between add_mutations_step() and zap_individuals_step()

void StepButtonCB(FL_OBJECT *p0, long p1) {
  static int step = 1;

  if (step % 2)
    add_mutations_step();
  else
    zap_individuals_step();

  step += 1;
}

