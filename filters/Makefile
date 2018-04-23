ifeq ($(PTHREADS),1)
	TYPE=-DWITH_PTHREADS
	CC 		= gcc
	CFLAGS 	+= $(TYPE) -g -Wall -Wextra
	LFLAGS	= -lpthread -lm -lrt
endif
ifeq ($(MPI),1)
	TYPE=-DWITH_MPI
	CC 		= mpicc
	CFLAGS 	+= $(TYPE) -g -Nmpi -Wall -Wextra
	LFLAGS	= -lpthread -lrt -lm
endif

all: blurc thresc

clean:
	-$(RM) *.o blurc thresc *.ppm *-results.png

blurc: ppmio.o gaussw.o blurfilter.o blurmain.o image_utils.o
	${CC} -o $@ $^ $(LFLAGS)

thresc: thresmain.o ppmio.o thresfilter.o
	${CC} -o $@ thresmain.o ppmio.o thresfilter.o $(LFLAGS)

arc:
	tar cvf filters.tar.gz *.c *.h *.f90 Makefile Makefile.sgi
