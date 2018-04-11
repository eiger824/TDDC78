CC 		= mpicc
CFLAGS 	+= $(TYPE) -g
LFLAGS	= -lpthread -lrt -lm

TYPE=-DWITH_MPI

ifeq ($(MPI),0)
	TYPE=-DWITH_OPENMP
endif

all: blurc thresc

clean:
	-$(RM) *.o blurc thresc

blurc: ppmio.o gaussw.o blurfilter.o blurmain.o
	${CC} -o $@ ppmio.o gaussw.o blurfilter.o blurmain.o $(LFLAGS)

thresc: thresmain.o ppmio.o thresfilter.o
	${CC} -o $@ thresmain.o ppmio.o thresfilter.o $(LFLAGS)

arc:
	tar cf - *.c *.h *.f90 Makefile Makefile.sgi|gzip - > filters.tar.gz
