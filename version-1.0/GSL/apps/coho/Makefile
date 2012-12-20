
# Coho salmon mutation simulation						
#                                                                            	
# The source files are								
#										
#  mas.C	the main simulation program					
#  scmas.C	sequential simulator source					
#  coho.*	the definition of class Coho, derived from the Population	
#		  class of the genetic simulation library			
#										

GSL=../..
GSLINCL=$(GSL)/include
GSLLIB=$(GSL)/lib

# C++ compiler flags for Sun:
CCFLAGS=-O3 -I$(GSLINCL)

# libraries:
SUNLIBS=-lgen -lm

# linker flags:
SUNLDFLAGS=-L$(GSLLIB)

scmas:	scmas.o coho.o
	CC $(SUNLDFLAGS) -o scmas scmas.o coho.o $(SUNLIBS)

scmas.o:	scmas.C coho.h cparams.h
coho.o:	coho.C coho.h

clean:
	/bin/rm -f core *.o *~ *.bak

