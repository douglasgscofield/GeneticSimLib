
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

//
// individual.h -- representation of a simulated individual
//

#ifndef _individual_h_
#define _individual_h_

#include <iostream>

class Genome;
class GBase;

enum sex_t {NONE, FEMALE, MALE};

class IndividualParamBlock {
};

class Individual {
public:
  Individual();
  ~Individual();
  static void reset_class();
  static void set_parameters(IndividualParamBlock *pb);
  void print(std::ostream &s);
  void set_sex(sex_t sex);
  sex_t get_sex();
  int get_id();
  Genome *const genes;
  void copy_genes(GBase *gb, int n);
protected:
  sex_t d_sex;
  int d_id;
  static int id;
};

inline std::ostream &operator <<(std::ostream &sout, Individual &x) {
  x.print(sout);
  return sout;
}

// Class ISet -- set of individuals.  

struct ILink {
  ILink *next;
  Individual *p;
};

class ISet {
public:
  ISet();
  ~ISet();
  void insert(Individual *x);
  Individual *remove(int i);
  Individual *select();
  Individual *operator [](int i);
  int size();
private:
  ILink *d_first;
  int d_size;
  static int count;
};


#endif  // _individual_h_
