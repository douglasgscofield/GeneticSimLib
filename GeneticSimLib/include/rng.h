
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/

/*                                                                           
  Random number generators.                                                  
                                                                             
  John Conery                                                                
  University of Oregon                                                       
                                                                             
  There are two ways to obtain a random deviate using procedures define      
  in this module.                                                            
                                                                             
  The first is to call a stand-alone function and pass it all the            
  required parameters.  For example, to obtain a number from a normal        
  distribution with mean M and standard deviation S, call rnormal(M,S).      
                                                                             
  The second technique is to create a random number generator object         
  with the parameters you want, and then use an operator to invoke a         
  method of the class to obtain a random deviate.  To obtain a normally      
  distributed numbrt with this technique, define a NormalRNG object, e.g.    
                                                                             
      NormalRNG r(M,S);                                                      
                                                                             
  and then every invocation of ++r will return a new number.                 
                                                                             
  This technique has two advantages over calling the function call:          
  complicated initialization code is executed just once,                     
  by the class constructor, and there is common interface for all            
  real-valued random deviates, so they can be implemented by virtual         
  operators.                                                                 
                                                                             
  The syntax used to communicate with an RNG is based on a standard          
  syntax for interacting with "iterators".  From this point of view,         
  a random distribution is a collection of values, and the object            
  iterates over that collection to return values from the distribution.      
  The iterator is initialized with parameters of the distribution,           
  e.g the mean and standard deviation of a Noraml distribution,              
  and then the iterator will cycle infinitely over the entire sequence.      
                                                                             
  The current distributions and the values that are either passed in a       
  stand-alone function call or used to initialize an RNG object are:         
                                                                             
    BinomialRNG(p,n) -- an integer between 0 and n, corresponding            
      to the number out of n trials that succeed, where each trial           
      has probability p.                                                     
                                                                             
    ExponentialRNG() -- a positive real value with exponentially             
      decreasing probability of higher values.                               
                                                                             
    Gamma(n) -- waiting time to the nth event in a Poison process            
                                                                             
    LogNormalRNG(a,b) -- random values such that the logarithms are	     
      normally distributed; a is the mean and b the standard deviation	     
      of the actual distribution.  This RNG estimates the mean and	     
      standard deviation of the underlying normal; use a second form,	     
      LogNormalLogRNG, if the actual mean and standard deviation are known.  
									     
    NormalRNG(a,b) -- a real value from a Gaussian distribution with         
      mean m and standard deviation s                                        
                                                                             
    PoissonRNG(lambda) -- an integer with expected value lambda and          
      exponentially decreasing probability of higher values.                 
                                                                             
    UniformRNG(a,b) -- a real number in the range [a,b]                      
                                                                             
  NOTE: the integer valued distributions (Poisson and Binomial) return       
  values of type double, just like the other distributions, but the          
  values returned are always whole numbers (e.g. 2.0, 5.0, ...)              
                                                                             
  In addition to the real-valued distributions listed above, there are       
  two stand-alone functions that return random integers:                     
                                                                             
    unsigned int rword() -- an integer in the range 0..2^N-1, where N is     
      the number of bits in an unsigned int (probably 32).                   
                                                                             
    unsigned int rword(n) -- an integer in the range 0..n-1 (simply the      
      value of (rword() % n)                                                 
                                                                             
  The underlying random sequences are produced by the "rand48" package of    
  mixed congruential generators.  By default all objects generated by any    
  concrete class derived from RNG will share the same underlying rand48      
  sequence, which will be automatically initialized from the system clock.   
  It is possible to call a base class method to define a particular seed     
  or to give each RNG object its own sequence.                               
                                                                             
  The definitions of the Binomial, Exponential, Normal, and Poisson          
  generators are all from Numerical Recipes in C by Press, et al.  In the    
  implementations below, I made no attempt to optimize any equations on the  
  right hand sides of assignment statements, and for the most part the       
  variable names are the same as those used by Press et al.  I did move      
  initialization code to constructor functions so it is executed only once,  
  and not once per call, and in some cases constant calculations have been   
  moved outside of loops.                                                    
                                                                             
  The LogNormal and Gamma generators are based on the descriptions from      
  "The Art of Computer Systems Performance Analysis" by Raj Jain             
  (J. Wiley and Sons, 1991).                                                 
                                                                             
  The syntax for iterators (using the ++ operator to obtain the next value   
  from a sequence) is from "Large Scale C++ Software Design" by J. Lakos     
  (Addison-Wesley, 1996).  The idea of a random number generator object and  
  the use of iterators to obtain values from the objects is original to this 
  implementation.                                                            
                                                                             
*/

#ifndef _rng_h_
#define _rng_h_

#include <stddef.h>

// A seed object represents the current state of a random sequence.
// Since we use the rand48 set of system functions, the state is
// 48 bits stored in an array of 3 unsigned shorts.  Ideally the
// state would be private, but that would hide the state from
// the random functions, too....

class RNGSeed {
public:
  RNGSeed();			// creates seed val from system clock
  RNGSeed(unsigned int x0);	// creates seed from integer
  virtual ~RNGSeed();
  unsigned short *x;	   	// 48 bits of state for erand48 et al
};

// Stand-alone functions.  The last parameter in each of these
// functions has a default value that points to a common seed
// defined in the implementation file, so the third parameter
// is optional in all these functions.

double rbinomial(double p, int n, RNGSeed *s = NULL);
double rexponential(RNGSeed *s = NULL);
double rgamma(double a, double b, RNGSeed *s = NULL);
double rlognormal(double mean, double stddev, RNGSeed *s = NULL);
double rlognormallog(double mean, double stddev, RNGSeed *s = NULL);
double rnormal(double mean, double stddev, RNGSeed *s = NULL);
double rpoisson(double lambda, RNGSeed *s = NULL);
double runiform(double lower, double upper, RNGSeed *s = NULL);
unsigned int rword(RNGSeed *s = NULL);
unsigned int rword(int n, RNGSeed *s = NULL);

// Random number generator base class.  RNG is derived from RNGSeed
// in order to give the random functions direct access to the state 
// of the seed.

class RNG {
public:
  virtual ~RNG();
  virtual double operator ++() =0;
  operator const double();
  RNGSeed *get_seed();			// return pointer to seed
  RNGSeed *save_seed();			// copy current state of seed
  void set_seed(RNGSeed *px);		// replace current seed with px
  void restore_seed(RNGSeed *ps);	// copy state back
  void inc_seed(int n);			// add n to the seed
protected:
  RNGSeed *seed;			// seed state for this object
  double prev;				// last value generated by this RNG
};

// Derived classes -- the last parameter of each constructor will
// default to a shared seed state defined in the implementation file.
// The other parameters define the distribution, e.g. the mean and
// standard deviation of a normal distribution.

class BinomialRNG : public RNG {
public:
  BinomialRNG(double p, int n, RNGSeed *s = NULL);
  double operator ++();
private:
  int n;
  double p;
  double xp;
  double b;
  double m;
  double g;
  double logp;
  double logq;
  double sq;
};

typedef enum {CDF_OK, CDF_OPEN_ERR, CDF_FORMAT_ERR} CDFStatus;

class CDFRNG : public RNG {
public:
  CDFRNG(const char* CDFFileName, RNGSeed *s = NULL);
  ~CDFRNG();
  int bad();
  double min();
  double max();
  double operator ++();
private:
  double *xv;
  int *cdf;
  int N;
  int rmax;
  CDFStatus err;
};

class ExponentialRNG : public RNG {
public:
  ExponentialRNG(RNGSeed *s = NULL);
  double operator ++();
};

class GammaRNG : public RNG {
public:
  GammaRNG(double a, double b, RNGSeed *s = NULL);
  double operator ++();
private:
  double a;
  double b;
};

class NormalRNG : public RNG {
public:
  NormalRNG(double mean, double stddev, RNGSeed *s = NULL);
  double operator ++();
private:
  double mean;
  double stddev;
};

class LogNormalRNG : public RNG {
public:
  LogNormalRNG(double mean, double stddev, RNGSeed *s = NULL);
  ~LogNormalRNG();
  double operator ++();
private:
  NormalRNG *r;
};

class LogNormalLogRNG : public RNG {
public:
  LogNormalLogRNG(double mean, double stddev, RNGSeed *s = NULL);
  ~LogNormalLogRNG();
  double operator ++();
private:
  NormalRNG *r;
};

class PoissonRNG : public RNG {
public:
  PoissonRNG(double lambda, RNGSeed *s = NULL);
  double operator ++();
private:
  double lambda;
  double sq;
  double alxm;
  double g;
};

class UniformRNG : public RNG {
public:
  UniformRNG(double lower, double upper, RNGSeed *s = NULL);
  double operator ++();
private:
  double lower;
  double range;
};

#endif  // _rng_h_
