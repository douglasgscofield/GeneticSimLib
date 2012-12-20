
/*--------------------------------------------------------------*/
/*  Copyright (c) 1997 by the University of Oregon.		*/
/*  See the COPYRIGHT file in this directory for permission	*/
/*  to use and distribute this software.			*/
/*--------------------------------------------------------------*/


/*
  Demo program for random number generator class
*/

#include "rng.h"

#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdlib.h>

#include "forms.h"
#include "rgui.h"

FD_rgui *the_gui;

#define NB 100
#define PROGRESS_THRESHOLD 100000

#define CHOICES "Binomial | CDF | Exponential | Gamma | LogNormal | LogNormalLog | Normal | Poisson | Uniform"
typedef enum {Binomial=1, CDF, Exponential, Gamma, LogNormal, LogNormalLog, Normal, Poisson, Uniform} choice;

/*
  When user clicks the "plot" button, read the values in the distribution
  choice and the A and B entries to create a random number generator.
*/

RNG *CreateRNG() {
  FL_OBJECT *ap = the_gui->AEntry;
  FL_OBJECT *bp = the_gui->BEntry;
  FL_OBJECT *cp = the_gui->DistributionChoice;
  char buf[80];
  char *cdffn;
  double a;
  double b;
  int i;
  RNG *r;

  switch (fl_get_choice(cp)) {
  case Uniform:
    a = atof(fl_get_input(ap));
    b = atof(fl_get_input(bp));
    r = new UniformRNG(a,b);
    break;
  case Binomial:
    a = atof(fl_get_input(ap));
    i = atoi(fl_get_input(bp));
    r = new BinomialRNG(a,i);
    break;
  case CDF:
    CDFRNG *cdfr;
    cdffn = (char *)fl_show_input("PDF File Name","wdist.dat");
    if (cdffn == NULL)
      return NULL;
    cdfr = new CDFRNG(cdffn);
    if (cdfr->bad()) {
      fl_show_message("Can't open PDF file",cdffn,"");
      return NULL;
    }
    sprintf(buf,"%f",cdfr->min());
    fl_set_input(ap,buf);
    sprintf(buf,"%f",cdfr->max());
    fl_set_input(bp,buf);
    r = cdfr;
    break;
  case Exponential:	
    r = new ExponentialRNG();
    break;
  case Gamma:
    a = atof(fl_get_input(ap));
    b = atof(fl_get_input(bp));
    r = new GammaRNG(a,b);
    break;
  case LogNormal:
    a = atof(fl_get_input(ap));
    b = atof(fl_get_input(bp));
    r = new LogNormalRNG(a,b);
    break;
  case LogNormalLog:
    a = atof(fl_get_input(ap));
    b = atof(fl_get_input(bp));
    r = new LogNormalLogRNG(a,b);
    break;
  case Normal:
    a = atof(fl_get_input(ap));
    b = atof(fl_get_input(bp));
    r = new NormalRNG(a,b);
    break;
  case Poisson:
    a = atof(fl_get_input(ap));
    r = new PoissonRNG(a);
    break;
  }
  return r;
}

// This routine is called when the user clicks the plot button.  Make
// an RNG, sample it a bunch of times, plot the results.

void PlotButtonCB(FL_OBJECT *p0, long p1) {
  FL_OBJECT *np = the_gui->NEntry;
  FL_OBJECT *sp = the_gui->ProgressSlider;
  FL_OBJECT *cp = the_gui->DistributionChoice;
  FL_OBJECT *xyp = the_gui->xyplot;
  RNG *rng = CreateRNG();
  int n = atoi(fl_get_input(np));
  double r;
  int i, j;
  int distrib[NB];
  float xv[NB];
  float yv[NB];
  int di;
  float dx;
  unsigned long xl;
  double sum, sumsq;
  static char *mstring = new char[80];
  static char *sdstring = new char[80];
  int ymax;
  double xmin, xmax;

  if (rng == NULL)
    return;

  for (i=0; i<NB; i++)
    distrib[i] = 0;
  sum = sumsq = 0.0;
  ymax = 0;
  
  if (n >= PROGRESS_THRESHOLD) {
    di = n/100;			// update slider every 'di' points
    fl_set_slider_value(sp,0.0);
    fl_show_object(sp);
  }
  else
    di = 0;			// n too small; don't draw slider

  xmin = atof(fl_get_input(the_gui->MinEntry));
  xmax = atof(fl_get_input(the_gui->MaxEntry));

  for (i=0; i<n; i++) {
    r = ++(*rng);		// get a sample from current distribution
    sum += r;
    sumsq += r*r;

    j = fl_get_choice(cp);
    if ((j == Binomial) || (j == Poisson))
      r *= 10.0;		// spread out integer values so we can see them
    
    if (j == CDF)		// 'dx' is a value added to each sample
      dx = 0.0005;		// the CDF demo RNG returns values that are exactly
    else			// on the divide between two bins, and roundoff
      dx = 0.0;			// errors lead to funny looking plots....

    j = floor((r+dx-xmin)*NB/(xmax-xmin));
    
    if ((j >= 0) && (j < NB)) {
      distrib[j] += 1;
      if (distrib[j] > ymax)
	ymax = distrib[j];
    }
    if (di && ((i%di)==0))
      fl_set_slider_value(sp,(double)i/(double)n);
  }

  fl_hide_object(sp);		// all done; put away the slider

  // Bins filled; transfer to x and y vectors and plot.  The xyplot
  // widget will scale x and y automatically.

  fl_delete_xyplot_text(xyp,mstring);
  fl_delete_xyplot_text(xyp,sdstring);

  for (i=0; i<NB; i++) {
    xv[i] = (float)i*(xmax-xmin)/NB + xmin;
    yv[i] = (float)distrib[i];
  }

  if (xmin == xmax)
    fl_set_xyplot_xbounds(xyp,xmin+1.0,xmax); 	// auto-scale if min>max
  else
    fl_set_xyplot_xbounds(xyp,xmin,xmax);
  fl_set_xyplot_data(xyp, xv, yv, NB, "", "", "");

  // compute, display sample mean and sample standard deviation

  double dn = double(n);
  double mean = sum/dn;
  double sd = sqrt(dn*((sumsq/dn)-mean*mean)/(dn-1));

  sprintf(mstring,"Sample Mean =  %f",mean);
  sprintf(sdstring,"Standard Dev = %f",sd);

  // The (x,y) coords of the string are defined in terms of scaled plot
  // units.  If the right edge of the x scale is at 100, put the text
  // at 80.  Put the top string near the top of the plot, and the bottom
  // string 10% lower.

  double y1 = double(ymax) - .05*(double(ymax));
  double y2 = double(ymax) - .10*(double(ymax));

  fl_add_xyplot_text(xyp,(xmax-xmin)*0.8+xmin,y1,mstring,FL_ALIGN_CENTER,FL_BLACK);
  fl_add_xyplot_text(xyp,(xmax-xmin)*0.8+xmin,y2,sdstring,FL_ALIGN_CENTER,FL_BLACK);
  
  delete rng;			// some RNGs allocate internal tables....
}

/*
  Initialize attributes of the XYPlot
*/

void InitPlot() {
  FL_OBJECT *xyp = the_gui->xyplot;
  FL_OBJECT *np = the_gui->NEntry;
  FL_OBJECT *sp = the_gui->ProgressSlider;
  FL_OBJECT *minp = the_gui->MinEntry;
  FL_OBJECT *maxp = the_gui->MaxEntry;

  fl_set_xyplot_fontsize(xyp,FL_NORMAL_SIZE);
  fl_set_xyplot_fontstyle(xyp,FL_FIXED_STYLE);
  fl_set_input(np,"10000");
  fl_set_input(minp,"0");
  fl_set_input(maxp,"100");
  fl_hide_object(sp);
}

// The user has changed the distribution.  Figure out which one they
// want to use now, update the labels on the arguments, and install
// default arguments for that type of distribution.

void DistributionChoiceCB(FL_OBJECT *p0, long p1) {
  FL_OBJECT *ap = the_gui->AEntry;
  FL_OBJECT *bp = the_gui->BEntry;

  switch (fl_get_choice(p0)) {
  case 0:
    cout << "No choice" << endl;
    break;
  case Uniform:			// uniform in the interval [a,b]
    fl_set_object_label(ap,"Minimum");
    fl_set_input(ap,"10.0");
    fl_set_object_label(bp,"Maximum");
    fl_set_input(bp,"90.0");
    fl_show_object(ap);
    fl_show_object(bp);
    break;
  case Binomial:		// binomial
    fl_set_object_label(ap,"p");
    fl_set_input(ap,"0.5");
    fl_set_object_label(bp,"n");
    fl_set_input(bp,"10.0");
    fl_show_object(ap);
    fl_show_object(bp);
    break;
  case CDF:			// CDF -- file will set range in a, b entries
    fl_set_object_label(ap,"lo");
    fl_set_input(ap,"<>");
    fl_set_object_label(bp,"hi");
    fl_set_input(bp,"<>");
    fl_show_object(ap);
    fl_show_object(bp);
    break;
  case Exponential:		// exponential; no params
    fl_hide_object(ap);
    fl_hide_object(bp);
    break;
  case Gamma:			// Time to wait for N events
    fl_set_object_label(ap,"Mean");
    fl_set_input(ap,"1");
    fl_set_object_label(bp,"Std Dev");
    fl_set_input(bp,"1");
    fl_show_object(ap);
    fl_show_object(bp);
    break;
  case LogNormal:		// lognormal with mean = a and std dev = b
    fl_set_object_label(ap,"Mean");
    fl_set_input(ap,"50.0");
    fl_set_object_label(bp,"Std Dev");
    fl_set_input(bp,"10.0");
    fl_show_object(ap);
    fl_show_object(bp);
    break;
  case LogNormalLog:		// lognormal with mean = exp(a) and std dev = exp(b)
    fl_set_object_label(ap,"Log Mean");
    fl_set_input(ap,"3.0");
    fl_set_object_label(bp,"Log Std Dev");
    fl_set_input(bp,"0.5");
    fl_show_object(ap);
    fl_show_object(bp);
    break;
  case Normal:			// normal with mean = a and std dev = b
    fl_set_object_label(ap,"Mean");
    fl_set_input(ap,"50.0");
    fl_set_object_label(bp,"Std Dev");
    fl_set_input(bp,"10.0");
    fl_show_object(ap);
    fl_show_object(bp);
    break;
  case Poisson:			// poisson with lambda = a
    fl_set_object_label(ap,"Lambda");
    fl_set_input(ap,"1.0");
    fl_show_object(ap);
    fl_hide_object(bp);
    break;
  }
}

/*
  Set the choice box text, sets the font, etc.
*/

void InitDistributionChoice() {
  FL_OBJECT *cp = the_gui->DistributionChoice;

  fl_clear_choice(cp);
  fl_addto_choice(cp,CHOICES);
  fl_set_choice_fontsize(cp,FL_NORMAL_SIZE);
  fl_set_choice_fontstyle(cp,FL_BOLD_STYLE);
  fl_set_choice(cp,Normal);
  DistributionChoiceCB(cp,0);
}

main(int argc, char *argv[]) {
  fl_initialize(&argc,argv,"Randoms",0,0);
  the_gui = create_form_rgui();
  InitDistributionChoice();
  InitPlot();
  fl_show_form(the_gui->rgui,FL_PLACE_CENTER,FL_FULLBORDER,"Random Distributions");

  while (1)
    fl_do_forms();
}

void QuitButtonCB(FL_OBJECT *p0, long p1) {
  exit(0);
}

