
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

//
// strand.C -- vector representation of strand
//

#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <iostream.h>

#include "strand.h"
#include "rng.h"


// Each vector entry is a single locus.  In this version, we're
// looking for the most compact representation, so we store only
// the locus index and the two mutation values.

// Since the strand is a vector of loci, don't define a constructor
// (otherwise it will be called for each vector element when a 
// strand is made....) [This could be a struct, but it's defined
// to be a class in the .h file]

class Locus {
public:
  int x;			// locus index (0..lmax-1)
  mutation_t s[2];		// the two mutations
};

int Strand::lmax = 10000;	// default number of loci per genotype
int Strand::lseg = 10;		// number of entries per vector 

// A strand is a vector of loci along with a count of the number
// of slots used and the total number of slots.

Strand::Strand() {
  v = NULL;
  nl = 0;
  ml = lseg;
}

Strand::~Strand() {
  if (v) {
    //    cerr << "deleting " << (void *)v << " [destructor]" << endl;
    delete v;
  }
}

// The copy and assignment operators are not allowed (declared private
// and not implemented):

Strand::Strand(const Strand &s) {
}

Strand Strand::operator=(const Strand &s) {
  return s;
}

void Strand::set_chromosome_length(int n) {
  lmax = n;	
  lseg = 10;			// ?? how can we estimate this from n?
}

// Add a mutation with effect sj at locus x.  0 <= y <= 1 indicates
// indicates which of the two strands gets the mutation.

fitness_t Strand::set_gene(int x, int y, mutation_t sj) {
  int i;
  int match;
  fitness_t wi, wj;
  mutation_t si;

  if (v == NULL) {
    v = new Locus[lseg];
    //    cerr << "new " << (void *)v << " [set_gene]" << endl;
  }

  for (i = 0, match = 0; i < nl; i++) {
    if (v[i].x == x) {
      match = 1;
      break;
    }
    if (x < v[i].x)
      break;
  }

  // If there is already a mutation at this locus (match == 1), 
  // see if the new mutation overwrites an existing one or if it is
  // for the other strand, then save the mutation.  If there is not
  // a mutation here yet, move the others down to make room for a
  // new locus in this strand and update the mutation count.

  if (match) {
    wi = locus_fitness(v[i].s[0],v[i].s[1]);
    if ((si = v[i].s[y]) != 0.0)
      v[i].s[y] = sj + si - sj*si;
    else
      v[i].s[y] = sj;
  }
  else {
    wi = 1.0;
    //    cerr << "set: ml = " << ml << ", nl = " << nl << ", i = " << i << endl;
    move_down(i);
    v[i].s[0] = v[i].s[1] = 0.0;
    v[i].x = x;
    v[i].s[y] = sj;
    nl += 1;
  }

  // compute new fitness and return relative decrease in fitness
  // from the new mutation at this location

  wj = locus_fitness(v[i].s[0],v[i].s[1]);
  //  cout << "set: w = " << wj/wi << endl;
  return wj/wi;
}

// We've run out of room in v.  Allocate a new, larger vector, copy
// the current loci, and toss out the old vector.

void Strand::extend_vector() {
  Locus *vt = v;

  v = new Locus[ml+lseg];
  //  cerr << "new " << (void *)v << " [extend]" << endl;
  //  cerr << "copying " << ml*sizeof(Locus)
  //       << " from " << (void *)vt << " to " << (void *)v << endl;
  memmove(v,vt,ml*sizeof(Locus));
  //  cerr << "deleting " << (void *)vt << " [extend]" << endl;
  delete vt;
  ml += lseg;
}

// Move current loci "down" one space to make room for a new one

void Strand::move_down(int i) {
  if (nl == ml)
    extend_vector();

  Locus *src = &v[i];
  Locus *dest = &v[i+1];

  int x = (nl-i);		// number of places to move
  //  cerr << "x = " << x
  //       << "#bytes = " << x*sizeof(Locus) 
  //       << " from " << (void *)src 
  //       << " to " << (void *)dest << endl;
  if (x) {
    memmove(dest,src,x*sizeof(Locus));
  }
}

// Return the value of mutation at locus x, strand y:

mutation_t Strand::get_gene(int x, int y) {
  int i;

  for (i = 0; i < nl; i++) {	// note: if v == NULL, nl will be 0
    if (v[i].x == x)
      return v[i].s[y];
    if (x < v[i].x)
      break;
  }

  return 0.0;
}

// Compute the fitness over all loci in this strand:

fitness_t Strand::strand_fitness() {
  fitness_t xw = 1.0;
  int i;

  for (i = 0; i < nl; i++)	// note: if v == NULL, nl will be 0
    xw *= locus_fitness(v[i].s[0],v[i].s[1]);
  
  return xw;
}

// Compute the fitness at locus x:

fitness_t Strand::get_fitness(int x) {
  int i;
  
  for (i = 0; i < nl; i++) {	// note: if v == NULL, nl will be 0
    if (v[i].x == x)
      return locus_fitness(v[i].s[0],v[i].s[1]);
    if (x < v[i].x)
      break;
  }
  
  return 1.0;			// no mutations here -- fitness = 1.0
}

// unravel(n,y) -- Create a list of loci from the loci in this object's
// strand.  Start with a random chromosome (top or bottom) and put the
// (non-zero) genes into chromosome 'y' of the new loci.  During the copy
// there will be 'n' cross-overs to the other input chromosome.

void sort(int *a, int n);	// simple insertion sort, defined below

Strand *Strand::unravel(int n, int y) {
  int i, j, k, x;
  static const int xmax = 100;	// max number of cross-overs per chromosome
  static int xx[xmax];		// array of cross-over points
  int cs;			// current strand (0 or 1);

  cs = rword(2);		// pick top or bottom strand as starting point

  for (i=0; i<n; i++)		// select n random cross-over points
    xx[i] = rword(lmax);	// (locus indices between 0 and max-1)
  xx[i] = lmax;
  sort(xx,n);			// cross-overs must be used in order

  // Create a return vector.  To be safe, make it as big as our own
  // vector:

  Locus *rv = new Locus[ml];
  //  cerr << "new " << (void *)rv << " [unravel]" << endl;

  // Scan the list of loci, starting from the head.  As we scan, apply the
  // cross-overs by switching strands at each cross-over index.  When the
  // gene on the current strand is non-zero, copy that locus to the return
  // vector

  //  cout << "-----" << endl;
  //  print(cout,10);
  //  cout << "unravel cs = " << cs << ", n = " << n << endl;

  i = 0;
  k = 0;
  for (j = 0; j < nl; j++) {
    while (v[j].x > xx[i]) {	// do cross-overs until we reach this locus
      cs ^= 1;			// (note: if A is 0 or 1, XOR(A,1) is the opposite)
      //      cout << "switch to " << cs << " " << xx[i] << " " << v[j].x << endl;
      i++;
    }
    if (v[j].s[cs] != 0.0) {	// is the non-zero mutation on current strand?
      rv[k].x = v[j].x;		// yes -- copy it
      rv[k].s[0] = 0.0;
      rv[k].s[1] = 0.0;
      rv[k].s[y] = v[j].s[cs];
      //      cout << "copy k = " << k << ", x = " << rv[k].x << ", s" << y << " = " << rv[k].s[y] << endl;
      k++;
    }
  }

  Strand *prs = new Strand;
  prs->v = rv;
  prs->nl = k;
  prs->ml = ml;

  return prs;
}

// merge(s1,s2) -- S1 and S2 are strands obtained by unraveling
// two parent strands.  Combine them into a single strand and return the
// fitness of the resulting vector.  

fitness_t Strand::merge(Strand *s1, Strand *s2) {
  fitness_t w = 1.0;		// return fitness
  int n1 = s1->nl;		// number of loci in strand
  int n2 = s2->nl;
  int i1 = 0;			// current locus indices
  int i2 = 0;
  int i = 0;
  Locus *v1 = s1->v;		// vectors from input strands
  Locus *v2 = s2->v;		
  int n = (s1->ml > s2->ml) ? s1->ml : s2->ml;

  // Allocate and initialize our own vector:

  if (v)
    delete v;
  v = new Locus[n];
  ml = n;
  nl = 0;

  // On each loop iteration, there are 5 cases:
  //  (a) no more loci in v1; copy the rest of v2 and exit
  //  (b) no more loci in v2; copy the rest of v1 and exit
  //  (c) the indices match; merge the two loci and continue
  //  (d) index in v1 < index in v2; copy from v1 and continue
  //  (e) index in v2 < index in v1; copy from v2 and continue

  //  cout << "s1: ";
  //  for (int j = 0; j < n1; j++)
  //    cout << " " << v1[j].x;
  //  cout << endl;
  //
  //  cout << "s2: ";
  //  for (j = 0; j < n2; j++)
  //    cout << " " << v2[j].x;
  //  cout << endl;

  fitness_t w0;

  while (1) {
    if (i1 >= n1) {
      while (i2 < n2) {
	//	cout << "[1] copy " << v2[i2].x;
	w0 = copy_locus(i++,v2,i2++);
	w *= w0;
	//	cout << " " << w0 << endl;
      }
      break;
    }
    if (i2 >= n2) {
      while (i1 < n1) {
	//	cout << "[2] copy " << v1[i1].x;
	w0 = copy_locus(i++,v1,i1++);
	w *= w0;
	//	cout << " " << w0 << endl;
      }
      break;
    }
    if (v1[i1].x == v2[i2].x) {
      //      cout << "[3] merge " << v1[i1].x;
      w0 = merge_locus(i++,v1,i1++,v2,i2++);
      w *= w0;
      //      cout << " " << w0 << endl;
      continue;
    }
    if (v1[i1].x < v2[i2].x) {
      //      cout << "[4] copy " << v1[i1].x;
      w0 = copy_locus(i++,v1,i1++);
      w *= w0;
      //      cout << " " << w0 << endl;
      continue;
    }
    //    cout << "[5] copy " << v2[i2].x;
    w0 = copy_locus(i++,v2,i2++);
    w *= w0;
    //    cout << " " << w0 << endl;
  }
  
  nl = i;
  delete s1;
  delete s2;

  //  cout << "** " << w << endl;
  return w;
}

fitness_t Strand::copy_locus(int i, Locus *p, int j) {
  if (i == ml)
    extend_vector();
  v[i].x = p[j].x;
  mutation_t s0 = v[i].s[0] = p[j].s[0];
  mutation_t s1 = v[i].s[1] = p[j].s[1];
  return locus_fitness(s0,s1);
}

// Since we know the input strands were created by 'unravel()',
// we know exactly one of the input vectors has a non-zero value...

fitness_t Strand::merge_locus(int i, Locus *p1, int j, Locus *p2, int k) {
  if (i == ml)
    extend_vector();
  v[i].x = p1[j].x;
  mutation_t s0 = v[i].s[0] = (p1[j].s[0] ? p1[j].s[0] : p2[k].s[0]);
  mutation_t s1 = v[i].s[1] = (p1[j].s[1] ? p1[j].s[1] : p2[k].s[1]);
  return locus_fitness(s0,s1);
}

// Print the loci in order; 'limit' is the maximum number to print

void Strand::print(ostream &sout, int limit) {
  int i;

  sout << "[" << nl << "/" << ml << "]";

  for (i = 0; (i < nl) && limit; i++, limit--) 
    sout << "  " << v[i].x << ":" << v[i].s[0] << "~" << v[i].s[1];
  
  if (i < nl)
    sout << "...";

  sout << endl;
}


// Strand iterator functions:

StrandIter::StrandIter(Strand &sp) {
  pv = sp.v;
  cur = 0;
  max = sp.nl;
}

int StrandIter::operator++() {
  if (cur < max) {
    cur += 1;
    return pv[cur-1].x;
  }
  else
    return Strand::lmax;
}

StrandIter::operator const int() const {
  if (cur < max)
    return pv[cur].x;
  else
    return Strand::lmax;
}

// local utility sort...

void sort(int *a, int n) {
  int i, j, t;

  for (i=1; i<n; i++) {
    t = a[i];
    for (j=i-1; j>=0 && (a[j] > t); j--)
      a[j+1] = a[j];
    a[j+1] = t;
  }
}
