
# Makefile for genetic simulation library

GSL=..
INCL=$(GSL)/include
LIB=$(GSL)/lib

#CPPFLAGS=-O3 -I$(INCL)
CPPFLAGS=-g -I$(INCL)

#CXXFLAGS=-O3 -m64 -I$(INCL)
CXXFLAGS=-g -m64 -I$(INCL)

SRCS=genome.cpp infgenome.cpp sparsegenome.cpp strand.cpp virtualgenome.cpp \
	individual.cpp generation.cpp population.cpp rng.cpp statistics.cpp gbase.cpp

OBJS=genome.o infgenome.o sparsegenome.o strand.o virtualgenome.o \
	individual.o generation.o population.o rng.o statistics.o gbase.o

libgen.a:	$(OBJS)
	ar -ru libgen.a $(OBJS)
	mv libgen.a $(LIB)

genome.o:	$(INCL)/genome.h genome.cpp
infgenome.o:	$(INCL)/genome.h $(INCL)/infgenome.h infgenome.cpp
sparsegenome.o:	$(INCL)/genome.h $(INCL)/sparsegenome.h $(INCL)/strand.h sparsegenome.cpp
virtualgenome.o:	$(INCL)/genome.h $(INCL)/virtualgenome.h virtualgenome.cpp
strand.o:	$(INCL)/strand.h strand.cpp
individual.o:	$(INCL)/individual.h individual.cpp
generation.o:	$(INCL)/generation.h generation.cpp
population.o:	$(INCL)/population.h population.cpp
rng.o:		$(INCL)/rng.h rng.cpp
statistics.o:	$(INCL)/statistics.h statistics.cpp
gbase.o:	$(INCL)/gbase.h gbase.cpp

clean:
	/bin/rm -f libgen.a $(LIB)/libgen.a core.* core *.o *~

