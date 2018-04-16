# TYPE=-DWITH_MPI

ifeq ($(PTHREADS),1)
	TYPE=-DWITH_PTHREADS
	CC 		= gcc
	CFLAGS 	+= $(TYPE) -g -Wall -Wextra
	LFLAGS	= -lpthread -lm -lrt
endif
ifeq ($(MPI),1)
	CC 		= mpicc
	CFLAGS 	+= $(TYPE) -g -Nmpi -Wall -Wextra
	LFLAGS	= -lpthread -lrt -lm
endif

all: blurc thresc

clean:
	-$(RM) *.o blurc thresc

blurc: ppmio.o gaussw.o blurfilter.o blurmain.o
	${CC} -o $@ ppmio.o gaussw.o blurfilter.o blurmain.o $(LFLAGS)

thresc: thresmain.o ppmio.o thresfilter.o
	${CC} -o $@ thresmain.o ppmio.o thresfilter.o $(LFLAGS)

arc:
	tar cvf filters.tar.gz *.c *.h *.f90 Makefile Makefile.sgi
