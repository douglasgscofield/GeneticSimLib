
#  mas -- mutation accumlation simulation

GSL=../..
GSLINCL=$(GSL)/include
GSLLIB=$(GSL)/lib

CPPFLAGS=-I$(GSLINCL) -O3
LDFLAGS=-L$(GSLLIB)
LIBS=-lm -lgen

mas:	mas.o
	g++ -o mas mas.o $(LDFLAGS) $(LIBS)

clean:
	/bin/rm -f mas core core.* *.o *~

