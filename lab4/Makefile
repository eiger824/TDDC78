#################### Variables ####################
CC 					= mpicc
CFLAGS 				= -c -Wall -Wextra -std=c99 -g
LDFLAGS 			= -lm

TARGET_SEQ			= sequential
TARGET_PARALLEL 	= parallel

OBJDIR				= obj
BINDIR 				= bin

OBJS_COMMON 		= ${OBJDIR}/physics.o \
					  ${OBJDIR}/log.o
OBJS				= ${OBJDIR}/main_sequential.o \
					  ${OBJDIR}/main_parallel.o \
					  ${OBJS_COMMON}
###################################################

all: ${OBJDIR} ${BINDIR} \
	${OBJS} \
	${BINDIR}/${TARGET_SEQ} \
	${BINDIR}/${TARGET_PARALLEL}

${OBJDIR}:
	test -d ${OBJDIR} || mkdir ${OBJDIR}

${BINDIR}:
	test -d ${BINDIR} || mkdir ${BINDIR}

${OBJDIR}/%.o: %.c
	${CC} ${CFLAGS} -o $@ $^

${BINDIR}/${TARGET_SEQ}: ${OBJDIR}/main_sequential.o ${OBJS_COMMON} 
	${CC} -o $@ $^ ${LDFLAGS} 

${BINDIR}/${TARGET_PARALLEL}: ${OBJDIR}/main_parallel.o ${OBJS_COMMON} 
	${CC} -o $@ $^ ${LDFLAGS} 

clean:
	rm -rf *.o *~ ${OBJDIR} ${BINDIR} *.swp
