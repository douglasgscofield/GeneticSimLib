
#  mas -- mutation accumlation simulation

GSL=../..
GSLINCL=$(GSL)/include
GSLLIB=$(GSL)/lib

CCFLAGS=-I$(GSLINCL) -O3
LDFLAGS=-L$(GSLLIB)
LIBS=-lm -lgen

mas:	mas.o
	CC -o mas mas.o $(LDFLAGS) $(LIBS)

clean:
	/bin/rm -f core *.o *~

