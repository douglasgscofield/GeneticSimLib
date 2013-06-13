
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/


//
// rng.cpp -- implementation of random number generator classes.
//
//  Note: the binomial and poisson functions use the math library
//  routine 'lgamma' to compute the Gamma function; according to the
//  man page it is not thread safe.  Change it to use lgamma_r and
//  recompile with _REENTRANT for threads based programs.
//

#include "rng.h"

#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <assert.h>
#include <fstream>
#include <string>
#include <cstring>

// Define a global seed to be used by stand-alone functions
// and objects that do not use a seed passed by the user.  Seeds
// such as this one that do not pass a value to the constructor
// are initialized using the system clock.

RNGSeed *default_seed = new RNGSeed;

// RNGSeed class implementation.  On 64-bit systems we can copy
// all 48 bits of state to the low 48 bits of the integer; on
// 32-bit systems the high order 16 bits will be lost, but 
// there will be no overflow errors (the bits are just shifted out).

// The default constructor initializes the seed with the system
// clock.  The other constructor assigns the seed to the user-
// supplied value.

// The 'inc' member function is useful in parallel programs where
// each thread or process gets the same time value from the system;
// to give each thread its own random sequence, have each thread
// add its thread ID to the seed.

RNGSeed::RNGSeed() {
  struct timeval t;
  gettimeofday(&t,NULL);
  x = new unsigned short[3];
  x[0] = 0;
  x[1] = (unsigned short)((t.tv_usec) >> 16);
  x[2] = (unsigned short)((t.tv_usec) & 0xFFFF);
}

RNGSeed::RNGSeed(unsigned int x0) {
  x = new unsigned short[3];
  // the (x0 >> 32) may cause a warning when unsigned int is 4 bytes
  x[0] = sizeof(x0) > 4 ? (x0 >> 32) & 0xFFFF : 0;
  x[1] = (x0 >> 16) & 0xFFFF;
  x[2] = x0 & 0xFFFF;
}

RNGSeed::~RNGSeed() {
  delete x;
}

// Base class operations...
//
// In all RNG classes, referring to a variable of type RNG returns the
// last value drawn from the distribution.
//

RNG::~RNG() {
}

RNGSeed *RNG::get_seed() {
  return seed;
}

RNGSeed *RNG::save_seed() {
  RNGSeed *px = new RNGSeed(0);
  px->x[0] = seed->x[0];
  px->x[1] = seed->x[1];
  px->x[2] = seed->x[2];
  return px;
}

void RNG::set_seed(RNGSeed *s) {
  seed = s;
}

void RNG::restore_seed(RNGSeed *s) {
  seed->x[0] = s->x[0];
  seed->x[1] = s->x[1];
  seed->x[2] = s->x[2];
}

void RNG::inc_seed(int n) {
  seed->x[2] += n;
}

RNG::operator const double() {
  return prev;
}

// The rword (random word) procedure.  Advances the sequence and
// returns an unsigned int constructed from the new seed state.
// For portability, only 32 bits of state is returned.....

unsigned int rword(RNGSeed *s) {
  unsigned short *xi = s ? s->x : default_seed->x;
  erand48(xi);
  return ((xi[1]<<16) + xi[2]);
}

unsigned int rword(int n, RNGSeed *s) {
  return rword(s) % n;
}

// --------- Derived Classes -------------

// All concrete classes derived from RNG need to provide a constructor
// and one operator.  The ++ operator advances the random sequence and
// returns the next deviate from the distribution.


// Binomial(p,n) -- the number of trials out of n that succeed, where
//    each trial has probability p.  This "helper function" does all
//    the work.

inline double binomial_f(int n, double p, double xp, double b, double m, double g,
			 double logp, double logq, double sq, unsigned short *xi) {

  int i;
  double y, em, t;

  do {
    do {
      y = tan(M_PI*erand48(xi));
      em = sq*y+m;
    } while (em < 0.0 || em >= (b+1.0));
    em = floor(em);
    t = 1.2*sq*(1.0+y*y)*exp(g-lgamma(em+1.0)-lgamma(b-em+1.0)+em*logp+(b-em)*logq);
  } while (erand48(xi) > t);

  return ((p == xp) ? em : (n - em));
}

// Binomial RNG constructor and operator

BinomialRNG::BinomialRNG(double pp, int pn, RNGSeed *s) {
  n = pn;
  p = (pp < 0.5 ? pp : (1.0 - pp));
  xp = pp;
  b = double(pn);
  m = b * p;
  g = lgamma(b+1.0);
  logp = log(p);
  logq = log(1.0-p);
  sq = sqrt(2.0*m*(1.0-p));
  seed = s ? s : default_seed;
}

double BinomialRNG::operator ++() {
  return prev = binomial_f(n,p,xp,b,m,g,logp,logq,sq,seed->x);
}

double rbinomial(double pp, int n, RNGSeed *s) {
  double p = (pp < 0.5 ? pp : (1.0 - pp));
  double b = double(n);
  double m = b * p;
  double g = lgamma(b+1.0);
  double logp = log(p);
  double logq = log(1.0-p);
  double sq = sqrt(2.0*m*(1.0-p));
  RNGSeed *sp = s ? s : default_seed;

  return binomial_f(n,p,pp,b,m,g,logp,logq,sq,sp->x);
}

// CDF(const std::string& f) -- read a probability density function (in the form
// of a histogram) from file f; use it to build a cumulative density
// function, and return random deviates from the CDF.


// Constructor function -- open a file, read the probability density
// function (a set of x-y pairs), and use it to create the cumulative
// density function.  Set 'err' to a nonzero value if the file can't
// be opened, if there is a format error on any line, or if the file
// does not contain the promised number of pairs.  The number of pairs
// is given on a comment line, e.g. "# N 100", that must occur before
// any data line.

//CDFRNG::CDFRNG(const std::string& CDFFileName, RNGSeed *s) {
//CDFRNG::CDFRNG(const char* CDFFileName, RNGSeed *s = NULL) {
CDFRNG::CDFRNG(const char* CDFFileName, RNGSeed *s) {
  std::ifstream in(CDFFileName);
  char line[100];
  int i;
  int n;
  char *nptr;
  float x;
  int y;

  if (in.bad()) {
    err = CDF_OPEN_ERR;
    return;
  }

  xv = NULL;
  cdf = NULL;

  while (in.getline(line,sizeof(line))) {
    if (line[0] == '\0')
      continue;
    if (line[0] == '#') {
      if ((nptr = strchr(line,'N')) && sscanf(nptr+1,"%d",&n)) {
	xv = new double[n+1];
	cdf = new int[n+1];
	i = 1;
	cdf[0] = 0;
	xv[0] = 0.0;
      }
      continue;
    }
    if ((i > n+1) || (cdf == NULL) || (xv == NULL)) {
      err = CDF_FORMAT_ERR;
      return;
    }
    if (sscanf(line,"%f %d",&x,&y) == 2) {
      xv[i] = double(x);
      cdf[i] = cdf[i-1] + y;
      i += 1;
    }
    else {
      err = CDF_FORMAT_ERR;
      return;
    }
  }
  if (i != n+1) {
    err = CDF_FORMAT_ERR;
    return;
  }

  N = n;
  rmax = cdf[n];
  err = CDF_OK;
  seed = s ? s : default_seed;
}

CDFRNG::~CDFRNG() {
  if (xv)
    delete xv;
  if (cdf)
    delete cdf;
}

int CDFRNG::bad() {
  return (err != CDF_OK);
}

double CDFRNG::min() {
  if (xv)
    return xv[1];
  else
    return 0.0;
}

double CDFRNG::max() {
  if (xv)
    return xv[N];
  else
    return 0.0;
}

// To draw a deviate, get a uniformly distributed integer between 1 
// and rmax, then find the i such that cdf[i-1] < r <= cdf[i].  To
// handle the special case when cdf[x] = cdf[x+1] (i.e. the original
// pdf is 0 at x) return the lowest value of x.

double CDFRNG::operator ++() {
  int i = N/2;			// initial search location
  int d = N/4;			// search increment
  int r = rword(rmax,seed)+1;	// the random integer to find in cdf[]
  
  while ((cdf[i] < r) || (cdf[i-1] > r)) {
    if (cdf[i] < r)
      i += d;
    else
      i -= d;
    if (d > 1)			// didn't find it -- halve the interval and continue
      d >>= 1;
  }

  while (cdf[i-1] == r)		// move to left of "plateau"
    i -= 1;

  return xv[i];			// return corresponding x value
}

// Exponential() -- a positive real value with exponentially
// decreasing probability of higher values

double rexponential(RNGSeed *s) {
  double x;
  RNGSeed *sp = s ? s : default_seed;
 
  do
    x = erand48(sp->x);
  while (x == 0.0);

  return -log(x);
}

ExponentialRNG::ExponentialRNG(RNGSeed *s) {
  seed = s ? s : default_seed;
}

double ExponentialRNG::operator ++() {
  return prev = rexponential(seed);
}

// Gamma distribution with mean a and std deviation b

inline double gamma_f(double a, double b, unsigned short *xi) {
  int i, bi;
  double x, y, xb, yb, bf;

  bi = b;			// separate b into integer part bi and
  bf = b - bi;			// fraction part bf

  if (bi) {
    x = 1.0;
    for (i = 0; i < bi; i++)
      x *= erand48(xi);
    x = -a * log(x);		// x is deviate from gamma with integer b
  }
  else
    x = 0.0;

  if (bf) {
    do {			// draw y from a beta distribution with 
      xb = erand48(xi);		// params bf and 1-bf (both less than 1)
      yb = erand48(xi);
      xb = pow(xb,1/bf);
      yb = pow(yb,1/(1-bf));
    } while (xb+yb > 1);
    y = xb/(xb+yb);
    y = a*y*rexponential();
  }
  else
    y = 0.0;

  return x + y;
}

GammaRNG::GammaRNG(double xa, double xb, RNGSeed *s) {
  a = (xb*xb)/xa;
  b = (xa*xa)/(xb*xb);
  seed = s ? s : default_seed;
}

double GammaRNG::operator ++() {
  return prev = gamma_f(a,b,seed->x);
}


double rgamma(double a, double b, RNGSeed *s) {
  RNGSeed *sp = s ? s : default_seed;
  return gamma_f(((b*b)/a),((a*a)/(b*b)),sp->x);
}

// LogNormal(m,s) and LogNormalLog(m,s) -- distribution with mean m
// and standard deviation s, where the logarithms of deviates are
// from a Gaussian distribution.  The first form estimates the mean
// and standard deviation of the underlying Gaussian, and the 
// deviates generated have mean m and standard deviation s.  In the
// second form, the user supplies the mean and standard of the
// underlying Gaussian, and the deviates returned have a mean that
// is approximately exp(m).

LogNormalRNG::LogNormalRNG(double mean, double stddev, RNGSeed *s) {
  double c = stddev/mean;
  double logstddev = sqrt(log(1+c*c));
  double logmean = log(mean) - log(1+c*c)/2;

  r = new NormalRNG(logmean,logstddev,s);
}

LogNormalRNG::~LogNormalRNG() {
  delete r;
}

double LogNormalRNG::operator ++() {
  return exp(++(*r));
}

double rlognormal(double mean, double stddev, RNGSeed *s) {
  double c = stddev/mean;
  double logstddev = sqrt(log(1+c*c));
  double logmean = log(mean) - log(1+c*c)/2;

  return exp(rnormal(logmean,logstddev,s));
}


LogNormalLogRNG::LogNormalLogRNG(double mean, double stddev, RNGSeed *s) {
  r = new NormalRNG(mean,stddev,s);
}

LogNormalLogRNG::~LogNormalLogRNG() {
  delete r;
}

double LogNormalLogRNG::operator ++() {
  return exp(++(*r));
}

double rlognormallog(double mean, double stddev, RNGSeed *s) {
  return exp(rnormal(mean,stddev,s));
}


// Normal(m,s) -- a sample from a Gaussian distribution with mean m
// and standard deviation s

inline double normal_f(double mean, double stddev, unsigned short *xi) {
  static int second = 0;
  static double saved;
  double dev, fac, r, v1, v2;

  if (second) {
    second = 0;
    dev = saved;
  }
  else {
    do {
      v1 = 2.0*erand48(xi)-1.0;
      v2 = 2.0*erand48(xi)-1.0;
      r = v1*v1+v2*v2;
    } while (r >= 1.0 || r == 0.0);
    /* (v1,v2) is a point inside the unit circle */
    fac = sqrt(-2.0*log(r)/r);
    saved = v1*fac;
    second = 1;
    dev = v2*fac;
  }

  return mean + (stddev * dev);
}

NormalRNG::NormalRNG(double m, double sd, RNGSeed *s) {
  mean = m;
  stddev = sd;
  seed = s ? s : default_seed;
}

double NormalRNG::operator ++() {
  return prev = normal_f(mean,stddev,seed->x);
}

double rnormal(double m, double sd, RNGSeed *s) {
  RNGSeed *sp = s ? s : default_seed;
  return normal_f(m,sd,sp->x);
}


// Poisson(a) -- Poisson distribution with expected value a

inline double poisson_f(double lambda, double sq, 
			double alxm, double g, unsigned short *xi) {
  int em;
  if (lambda < 12.0) {
    em = -1;
    double t = 1.0;
    do {
      em += 1;
      t *= erand48(xi);
    } while (t > g);
  }
  else {
    double t, u, y;
    do {
      do {
	y = tan(M_PI*erand48(xi));
	u = sq*y + lambda;
      } while (u < 0.0);
      em = floor(u+0.5);
      t = 0.9 * (1.0 * y*y) * exp(em*alxm-lgamma(em+1)-g);
    } while (erand48(xi) > t);
  }

  return double(em);
}

PoissonRNG::PoissonRNG(double a, RNGSeed *s) {
  lambda = a;
  sq = sqrt(2.0*lambda);
  alxm = log(lambda);
  if (lambda < 12.0)
    g = exp(-lambda);
  else
    g = lambda*alxm - lgamma(lambda+1.0);
  seed = s ? s : default_seed;
}

double PoissonRNG::operator ++() {
  return prev = poisson_f(lambda,sq,alxm,g,seed->x);
}

double rpoisson(double a, RNGSeed *s) {
  double lambda = a;
  double sq = sqrt(2.0*lambda);
  double alxm = log(lambda);
  double g; 
  RNGSeed *sp = s ? s : default_seed;

  if (lambda < 12.0)
    g = exp(-lambda);
  else
    g = lambda*alxm - lgamma(lambda+1.0);

  return poisson_f(lambda,sq,alxm,g,sp->x);
}

// Uniform(a,b) -- a uniformly varying value between a and b

inline double uniform_f(double lower, double range, unsigned short *xi) {
  return range*erand48(xi) + lower;
}

UniformRNG::UniformRNG(double a, double b, RNGSeed *s) {
  lower = a;
  range = b-a;
  seed = s ? s : default_seed;
}

double UniformRNG::operator ++() {
  return prev = uniform_f(lower,range,seed->x);
}

double runiform(double a, double b, RNGSeed *s) {
  double lower = a;
  double range = b-a;
  RNGSeed *sp = s ? s : default_seed;

  return uniform_f(lower,range,sp->x);
}

