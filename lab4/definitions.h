#ifndef _definitions_h
#define _definitions_h

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#include "coordinate.h"
#include "physics.h"

#define PI                      3.141592653

#define MAX_NO_PARTICLES        15000  /* Maximum number of particles/processor */
#define INIT_NO_PARTICLES       500    /* Initial number of particles/processor */
#define MAX_INITIAL_VELOCITY    50

#define DEFAULT_TSTEP           1.00

#define BOX_HORIZ_SIZE          10000.0
#define BOX_VERT_SIZE           10000.0
#define WALL_LENGTH             (2.0*BOX_HORIZ_SIZE+2.0*BOX_VERT_SIZE)

#define PARTICLE_BUFFER_SIZE    MAX_NO_PARTICLES/5
#define COMM_BUFFER_SIZE        5*PARTICLE_BUFFER_SIZE

/* MPI notation */
#define ROOT                    0

/* Filename to store the output data */
#define FILENAME                "results.txt"

typedef struct particle
{
    pcord_t  pcord;
    int ptype;        /* Used to simulate mixing of gases */ 
} particle_t;

typedef unsigned int uint ;

typedef enum neighbor_location
{
    LEFT = 0,
    TOP,
    RIGHT,
    BOTTOM
} nloc_t;

#endif

