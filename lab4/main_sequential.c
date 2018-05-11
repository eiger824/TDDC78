#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>


#include "coordinate.h"
#include "definitions.h"
#include "physics.h"


//Feel free to change this program to facilitate parallelization.

float rand1()
{
    return (float)( rand()/(float) RAND_MAX );
}

void init_collisions(bool *collisions, unsigned int max)
{
    for(unsigned int i=0;i<max;++i)
        collisions[i]=0;
}


int main(int argc, char** argv)
{
    unsigned int time_stamp = 0, time_max;
    float pressure = 0;
    int nr_particles = 8000;

    // parse arguments
    if (argc != 2)
    {
        fprintf(stderr, "Usage: %s simulation_time\n", argv[0]);
        fprintf(stderr, "For example: %s 10\n", argv[0]);
        exit(1);
    }

    time_max = atoi(argv[1]);

    /* Initialize */
    // 1. set the walls
    cord_t wall;
    wall.y0 = wall.x0 = 0;
    wall.x1 = BOX_HORIZ_SIZE;
    wall.y1 = BOX_VERT_SIZE;


    // 2. allocate particle bufer and initialize the particles
    pcord_t * particles = (pcord_t*) malloc(nr_particles * sizeof(pcord_t) );
    bool * collisions = (bool *) malloc (nr_particles * sizeof(bool) );

    srand( time(NULL) + 1234 );

    float r, a;
    for (int i = 0; i < nr_particles; i++)
    {
        // initialize random position
        particles[i].x = wall.x0 + rand1()*BOX_HORIZ_SIZE;
        particles[i].y = wall.y0 + rand1()*BOX_VERT_SIZE;

        // initialize random velocity
        r = rand1()*MAX_INITIAL_VELOCITY;
        a = rand1()*2*PI;
        particles[i].vx = r*cos(a);
        particles[i].vy = r*sin(a);
    }


    unsigned int p, pp;

    /* Main loop */
    for (time_stamp = 0; time_stamp < time_max; time_stamp++) // for each time stamp
    {
        init_collisions(collisions, nr_particles);

        for (p = 0; p < nr_particles; p++) // for all particles
        {
            if (collisions[p]) continue;

            /* check for collisions */
            for (pp = p + 1; pp < nr_particles; pp++)
            {
                if (collisions[pp]) continue;

                float t = collide(&particles[p], &particles[pp]);
                if (t != -1) // collision
                {
                    collisions[p] = collisions[pp] = 1;
                    interact(&particles[p], &particles[pp], t);
                    break; // only check collision of two particles
                }
            }
        }

        // move particles that has not collided with another
        for (p = 0; p < nr_particles; p++)
        {
            if (!collisions[p])
            {
                feuler(&particles[p], 1);

                /* check for wall interaction and add the momentum */
                pressure += wall_collide(&particles[p], wall);
            }
        }
    }


    printf("Average pressure = %f\n", pressure / (WALL_LENGTH*time_max));

    free(particles);
    free(collisions);

    return 0;

}

