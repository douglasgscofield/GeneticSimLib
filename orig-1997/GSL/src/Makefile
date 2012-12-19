
# Makefile for genetic simulation library

GSL=..
INCL=$(GSL)/include
LIB=$(GSL)/lib

CCFLAGS=-O3 -I$(INCL)
#CCFLAGS=-g -I$(INCL)

CXXFLAGS=-O3 -64 -I$(INCL)
#CXXFLAGS=-g -64 -I$(INCL)

SRCS=genome.C infgenome.C sparsegenome.C strand.C virtualgenome.C \
	individual.C generation.C population.C rng.C statistics.C gbase.C

OBJS=genome.o infgenome.o sparsegenome.o strand.o virtualgenome.o \
	individual.o generation.o population.o rng.o statistics.o gbase.o

libgen.a:	$(OBJS)
	ar -ru libgen.a $(OBJS)
	mv libgen.a $(LIB)

genome.o:	$(INCL)/genome.h genome.C
infgenome.o:	$(INCL)/genome.h $(INCL)/infgenome.h infgenome.C
sparsegenome.o:	$(INCL)/genome.h $(INCL)/sparsegenome.h $(INCL)/strand.h sparsegenome.C
virtualgenome.o:	$(INCL)/genome.h $(INCL)/virtualgenome.h virtualgenome.C
strand.o:	$(INCL)/strand.h strand.C
individual.o:	$(INCL)/individual.h individual.C
generation.o:	$(INCL)/generation.h generation.C
population.o:	$(INCL)/population.h population.C
rng.o:		$(INCL)/rng.h rng.C
statistics.o:	$(INCL)/statistics.h statistics.C
gbase.o:	$(INCL)/gbase.h gbase.C

clean:
	/bin/rm -f core *.o *~

