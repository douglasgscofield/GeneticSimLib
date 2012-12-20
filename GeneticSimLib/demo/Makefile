
# Makefile for demos and exercisers for genetic simulation library

CPP=g++

GSL=..
GSLINCL=$(GSL)/include
GSLLIB=$(GSL)/lib

# If you are going to compile the demos that use the XFORMS user
# interface package, edit these macros to point to the FORMS
# directory and X11 libraries on your system:

#FORMSDIR=/local/apps/xforms/FORMS
FORMSDIR=/opt/local/lib
#XINCLUDE=/usr/openwin/share/include
XINCLUDE=/opt/local/include
#XLIB=/lib/X11
XLIB=/opt/local/lib

CPPFLAGS=-I$(FORMSDIR) -I$(XINCLUDE) -I$(GSLINCL) -O3
LDFLAGS=-L$(FORMSDIR) -L$(XLIB) -L$(GSLLIB)
#LIBS=-lm -lforms -lX -lgen
LIBS=-lm -lforms -lX11 -lgen

all:	gdemo gdriver idemo rnddemo

gdemo:	gdemo.o ggui.o $(GSLLIB)/libgen.a
	$(CPP) -o gdemo gdemo.o ggui.o $(LDFLAGS) $(LIBS)

gdriver:	gdriver.o $(GSLLIB)/libgen.a
	$(CPP) -o gdriver gdriver.o $(LDFLAGS) $(LIBS)

idemo:	idemo.o $(GSLLIB)/libgen.a
	$(CPP) -o idemo idemo.o $(LDFLAGS) $(LIBS)

rnddemo:	rnddemo.o rgui.o $(GSLLIB)/libgen.a
	$(CPP) -o rnddemo rnddemo.o rgui.o $(LDFLAGS) $(LIBS)

# compile C programs with C++ compiler:

ggui.o:	ggui.c
	$(CPP) $(CPPFLAGS) -c ggui.c

rgui.o:	rgui.c
	$(CPP) $(CPPFLAGS) -c rgui.c

clean:
	/bin/rm -f core *~ *.o


