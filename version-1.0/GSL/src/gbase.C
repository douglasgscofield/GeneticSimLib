
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

/*
  gbase.C -- genotype database class

  The constructor opens a text file created by 'mkgbase' and reads a
  subset of the genotypes.  Member functions select an item in the
  database and return the values and locations of mutations in the
  currently selected item.
*/

#include <fstream.h>

#include "gsl.h"		// genetic simulation library 
#include "gbase.h"		// genotype database structure definition

// The default constructor (no parameters) is used by parallel applications.
// One PE will use the regular constructor to read the database from a
// file; the others will make empty GBase objects.  A derived class in
// the parallel app will have a method to copy a GBase object to an
// empty object.

GBase::GBase() {
  params = NULL;
  sv = NULL;	
  gv = NULL;
  gb_size = 0;
  cgx = 0;
  clx = 0;
}

int read_header(ifstream &f, gb_param p[]);
gb_mutation *read_fitness(ifstream &f, gb_param p[]);
gb_genotype *read_genotypes(ifstream &f, int &n, gb_param p[], gb_mutation sv[]);

GBase::GBase(int n, char *gfn) {
  ifstream f(gfn);
  char bigline[1024];

  if (f.bad()) 
    cerr << "GBase: Can't open " << gfn << endl;
  else {
    //    f.setbuf(bigline,sizeof(bigline));
    params = new gb_param[GB_NPARAMS];
    if (read_header(f,params)) {
      sv = read_fitness(f,params);
      gv = read_genotypes(f,n,params,sv);
      gb_size = n;
      cgx = 0;
      clx = 0;
    }
  }
}

int read_header(ifstream &f, gb_param p[]) {
  char line[80];
  int sizes[GB_NSEGS];
  int i;

  f.getline(line,sizeof(line));
  if (strcmp(line,GBTAG) != 0) {
    cerr << "GBase: incorrect tag in genotype database file" << endl;
    return 0;
  }

  for (i = 0; i < GB_NSEGS; i++)
    f >> sizes[i];			// not used in object...

  for (i = 0; i < GB_NPARAMS; i++)
    f >> p[i];

  return 1;
}

gb_mutation *read_fitness(ifstream &f, gb_param p[]) {
  int gl = p[GB_GL];
  gb_mutation *sv = new gb_mutation[gl];

  for (int i = 0; i < gl; i++)
    f >> sv[i];

  return sv;
}

// Return the location in array a (which currently has n items) where 
// we can find key k or would insert key k.

int loc(int *a, int n, int k) {
  int u = n;
  int l = 0;
  int m;

  while ((u-l) > 1) {
    m = (l+u)/2;
    if (k == a[m])
      return m;
    if (k < a[m])
      u = m;
    else
      l = m;
  }

  if (k <= a[l])
    return l;
  else
    return u;
}

// Create an array of N integers, each in the range [0..M-1], sorted
// in ascending order with no duplicates.  This method, which does a
// binary search followed by an insertion, will take a very long time
// if N is close to M....

int *make_sample(int N, int M) {
  int *a = new int[N];
  int i, j, k, l;

  a[0] = rword(M);
  j = 1;

  while (j < N) {
    k = rword(M);
    i = loc(a,j,k);		// find spot for k in a[0]..a[j-1]
    if (a[i] == k)
      continue;			// it's a duplicate; skip it
    for (l = j; l > i; l--)
      a[l] = a[l-1];
    a[i] = k;
    j += 1;
  }

  return a;
}

void skip(ifstream &f, int n) {
  int i, j, nm;
  gb_locus l;

  for (i = 0; i < n; i++) {
    f >> dec >> nm;
    for (j = 0; j < nm; j++)
      f >> hex >> l.x;
  }
}

void read_one(ifstream &f, gb_genotype *g) {
  int i, nm;
  gb_locus l;

  f >> dec >> nm;
  g->n = nm;
  g->m = new gb_locus[nm];
  for (i = 0; i < nm; i++) {
    f >> hex >> l.x;
    g->m[i] = l;
  }
}

void read_all(ifstream &f, int n, gb_genotype g[]) {
  int i;
  gb_genotype x;

  for (i = 0; i < n; i++) {
    read_one(f,&x);
    g[i] = x;
  }
}

// Read n random genotypes from the file f, put them in g

void read_selected(ifstream &f, int n, int ng, gb_genotype g[]) {
  int i, j;
  gb_genotype x;
  int *a = make_sample(n,ng);		// n random indices

  j = 0;				// j is index in file,
  for (i = 0; i < n; i++) {		// i is index in array
    skip(f,a[i]-j);
    read_one(f,&x);
    g[i] = x;
    j = a[i]+1;
  }
}

gb_genotype *read_genotypes(ifstream &f, int &n, gb_param p[], gb_mutation sv[]) {
  int ng = p[GB_NG];
  int i, j;
  int nm;
  gb_locus l;
  gb_genotype *g;

  if (ng < n)
    cerr << "GBase: request greater than GBase size" << endl;

  if (ng <= n)
    n = ng;			// note: n is a reference parameter...

  g = new gb_genotype[ng];

  if (ng == n)
    read_all(f,ng,g);
  else
    read_selected(f,n,ng,g);

  return g;
}

GBase::~GBase() {
  int n = params[GB_NG];
  int i;

  for (i = 0; i < n; i++)
    delete[] gv[i].m;

  delete[] gv;
  delete[] sv;
  delete[] params;
}

int GBase::size() {
  return gb_size;		// NOT params[GB_NG] -- that's the number of
				// genotypes in the file, not the number in the object
}

void GBase::select(int n) {
  cgx = (n % gb_size);		// makes sure index is in range
  clx = 0;
}

int GBase::num_mutations() {
  return gv[cgx].n;
}


int GBase::next_locus() {
  return (++clx < gv[cgx].n);
}

int GBase::locus_index() {
  return gv[cgx].m[clx].s.i;
}

int GBase::locus_type() {
  return gv[cgx].m[clx].s.f;
}

mutation_t GBase::locus_val(int n) {
  int ltype = gv[cgx].m[clx].s.f;
  if (((ltype == 1) && (n == 1)) || ((ltype == 2) && (n == 0)))
    return 0.0;
  return mutation_t(sv[gv[cgx].m[clx].s.i]);
}

mutation_t GBase::ms() {
  return params[GB_S];
}

mutation_t GBase::sds() {
  return params[GB_SDS];
}

double GBase::u() {
  return params[GB_U];
}

int GBase::gl() {
  return params[GB_GL];
}



