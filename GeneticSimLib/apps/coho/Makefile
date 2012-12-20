
# Coho salmon mutation simulation						
#                                                                            	
# The source files are								
#										
#  mas.cpp	the main simulation program					
#  scmas.cpp	sequential simulator source					
#  coho.*	the definition of class Coho, derived from the Population	
#		  class of the genetic simulation library			
#										

GSL=../..
GSLINCL=$(GSL)/include
GSLLIB=$(GSL)/lib

# C++ compiler flags for Sun:
#CPPFLAGS=-O3 -I$(GSLINCL)
CPPFLAGS=-g -I$(GSLINCL)

# libraries:
SUNLIBS=-lgen -lm

# linker flags:
SUNLDFLAGS=-L$(GSLLIB)

scmas:	scmas.o coho.o
	g++ $(SUNLDFLAGS) -o scmas scmas.o coho.o $(SUNLIBS)

scmas.o:	scmas.cpp coho.h cparams.h
coho.o:	coho.cpp coho.h

clean:
	/bin/rm -f scmas core.* core *.o *~ *.bak

