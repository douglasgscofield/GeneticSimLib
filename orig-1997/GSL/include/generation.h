
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

//
// generation.h -- a Generation is a collection of individuals born
// in the same reproductive cycle.
//

#ifndef _generation_h_
#define _generation_h_

class Individual;
class ISet;

class GenerationParamBlock {
};

class Generation {
public:
  Generation();
  ~Generation();
  void insert(Individual *x);
  Individual *remove(int i);
  Individual *select();
  Individual *operator [](int i);
  int size();
  static void reset_class();
  static void set_parameters(GenerationParamBlock *pb);
private:
  ISet *is;
};

#endif

