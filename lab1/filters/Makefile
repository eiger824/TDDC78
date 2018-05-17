OBJDIR = obj
BINDIR = bin

# For tracing with ITAC
CFLAGS += -tcollect -I$(VT_ROOT)/include
LFLAGS += -L$(VT_LIB_DIR) -lVT $(VT_ADD_LIBS)

ifeq ($(PTHREADS),1)
	CC 		= gcc
	CFLAGS 	+= -g -Wall -Wextra
	LFLAGS	+= -lpthread -lrt -lm
	# Blurfilter
	BLURFILTER_OBJ=blurfilter-pthreads.o
	BLURFILTER_SRC=blurfilter-pthreads.c
	BLURMAIN_OBJ=blurmain-pthreads.o
	BLURMAIN_SRC=blurmain-pthreads.c
	# Thresfilter
	THRESFILTER_OBJ=thresfilter.o
	THRESFILTER_SRC=thresfilter.c
	THRESMAIN_OBJ=thresmain-pthreads.o
	THRESMAIN_SRC=thresmain-pthreads.c
endif
ifeq ($(MPI),1)
	CC 		= mpicc
	CFLAGS 	+= -g -Nmpi -Wall -Wextra
	LFLAGS	+= -lpthread -lrt -lm
	# Blurfilter
	BLURFILTER_OBJ=blurfilter-mpi.o
	BLURFILTER_SRC=blurfilter-mpi.c
	BLURMAIN_OBJ=blurmain-mpi.o
	BLURMAIN_SRC=blurmain-mpi.c
	# Thresfilter
	THRESFILTER_OBJ=thresfilter.o
	THRESFILTER_SRC=thresfilter.c
	THRESMAIN_OBJ=thresmain-mpi.o
	THRESMAIN_SRC=thresmain-mpi.c
endif

all: test ${OBJDIR} ${BINDIR} ${BINDIR} ${BINDIR}/blurc ${BINDIR}/thresc

clean:
	rm -rf *.o blurc thresc *.ppm *-results.png ${OBJDIR} ${BINDIR} *.swp *~ results*.txt

test:
	@test -n "$(MPI)$(PTHREADS)"

${OBJDIR}:
	test -d ${OBJDIR} || mkdir ${OBJDIR}

${BINDIR}:
	test -d ${BINDIR} || mkdir ${BINDIR}

${OBJDIR}/%.o: %.c
	${CC} -c $(CFLAGS) $< -o $@

${BINDIR}/blurc: ${OBJDIR}/ppmio.o \
	${OBJDIR}/gaussw.o \
	${OBJDIR}/${BLURFILTER_OBJ} \
	${OBJDIR}/${BLURMAIN_OBJ} \
	${OBJDIR}/image_utils.o
	${CC} $^ $(LFLAGS) -o $@

${BINDIR}/thresc: ${OBJDIR}/ppmio.o \
	${OBJDIR}/${THRESFILTER_OBJ} \
	${OBJDIR}/${THRESMAIN_OBJ}
	${CC} $^ $(LFLAGS) -o $@

arc:
	tar jcvf filters-labs1-2-sanpa993-haoli436.tar.bz2 *.c *.h Makefile README scripts *.pdf 
