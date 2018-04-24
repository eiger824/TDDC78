OBJDIR = obj
BINDIR = bin

ifeq ($(PTHREADS),1)
	CC 		= gcc
	CFLAGS 	+= -g -Wall -Wextra
	LFLAGS	= -lpthread -lm -lrt
	BLURFILTER_OBJ=blurfilter-pthreads.o
	BLURFILTER_SRC=blurfilter-pthreads.c
	BLURMAIN_OBJ=blurmain-pthreads.o
	BLURMAIN_SRC=blurmain-pthreads.c
endif
ifeq ($(MPI),1)
	CC 		= mpicc
	CFLAGS 	+= -g -Nmpi -Wall -Wextra
	LFLAGS	= -lpthread -lrt -lm
	BLURFILTER_OBJ=blurfilter-mpi.o
	BLURFILTER_SRC=blurfilter-mpi.c
	BLURMAIN_OBJ=blurmain-mpi.o
	BLURMAIN_SRC=blurmain-mpi.c
endif

all: test ${OBJDIR} ${BINDIR} ${BINDIR} ${BINDIR}/blurc ${BINDIR}/thresc

clean:
	rm -rf *.o blurc thresc *.ppm *-results.png ${OBJDIR} ${BINDIR}

test:
	@test -n "$(MPI)$(PTHREADS)"

${OBJDIR}:
	test -d ${OBJDIR} || mkdir ${OBJDIR}

${BINDIR}:
	test -d ${BINDIR} || mkdir ${BINDIR}

${OBJDIR}/%.o: %.c
	${CC} -c $(CFLAGS) $< -o ${OBJDIR}/$@

${BINDIR}/blurc: ppmio.o \
	gaussw.o \
	${BLURFILTER_OBJ} \
	${BLURMAIN_OBJ} \
	image_utils.o
	${CC} -o $@ $^ $(LFLAGS)

${BINDIR}/thresc: thresmain.o ppmio.o thresfilter.o
	${CC} -o $@ thresmain.o ppmio.o thresfilter.o $(LFLAGS)

arc:
	tar cvf filters.tar.gz *.c *.h *.f90 Makefile
