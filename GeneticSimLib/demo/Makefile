
# Makefile for demos and exercisers for genetic simulation library

GSL=..
GSLINCL=$(GSL)/include
GSLLIB=$(GSL)/lib

# If you are going to compile the demos that use the XFORMS user
# interface package, edit these macros to point to the FORMS
# directory and X11 libraries on your system:

FORMSDIR=/local/apps/xforms/FORMS
XINCLUDE=/usr/openwin/share/include
XLIB=/lib/X11

CCFLAGS=-I$(FORMSDIR) -I$(XINCLUDE) -I$(GSLINCL) -O3
LDFLAGS=-L$(FORMSDIR) -L$(XLIB) -L$(GSLLIB)
LIBS=-lm -lforms -lX -lgen

all:	gdemo gdriver idemo rnddemo

gdemo:	gdemo.o ggui.o $(GSLLIB)/libgen.a
	CC -o gdemo gdemo.o ggui.o $(LDFLAGS) $(LIBS)

gdriver:	gdriver.o $(GSLLIB)/libgen.a
	CC -o gdriver gdriver.o $(LDFLAGS) $(LIBS)

idemo:	idemo.o $(GSLLIB)/libgen.a
	CC -o idemo idemo.o $(LDFLAGS) $(LIBS)

rnddemo:	rnddemo.o rgui.o $(GSLLIB)/libgen.a
	CC -o rnddemo rnddemo.o rgui.o $(LDFLAGS) $(LIBS)

# compile C programs with C++ compiler:

ggui.o:	ggui.c
	CC $(CCFLAGS) -c ggui.c

rgui.o:	rgui.c
	CC $(CCFLAGS) -c rgui.c

clean:
	/bin/rm -f core *~ *.o


