
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

// Genome.C -- procedures and data structures common to all genome classes.

// Implementation of the "virtual constructor" -- The examplar is an object 
// created by one of the derived types. The user must call the set_parameters()
// procedure of one of the derived classes to set the parameters of that 
// class and also to create an examplar object.  After the exemplar is 
// created, a new Genome object is created by calling make_genome(), which
// simply forwards the request to the new_genome() function of the
// current derived class.

#include "genome.h"
#include <stddef.h>

Genome *Genome::exemplar = NULL;

Genome *make_genome() {
  return Genome::exemplar->new_genome();
}

