#ifndef __DOUBLYLINKEDLIST_H
#define __DOUBLYLINKEDLIST_H

#include <stdbool.h>

#include "definitions.h"

typedef struct collisions_particles
{
    bool        collision;
    pcord_t  *  particle;
} part_coll_t;


typedef struct dll_node_type
{
    part_coll_t * p;
    struct dll_node_type *next;
    struct dll_node_type *prev;

} dll_node_t;

typedef struct dll_type
{
    dll_node_t *first;
    dll_node_t *last;

    int count;

} dll_t;

dll_t *dll_init(void);
void dll_destroy(dll_t * list);

void dll_insert_after(dll_t * list, dll_node_t *node, part_coll_t * p);
void dll_insert_before(dll_t * list, dll_node_t *node, part_coll_t * p);

void dll_insert_beginning(dll_t * list, part_coll_t * p); 
void dll_prepend(dll_t * list, part_coll_t * p);

void dll_insert_end(dll_t * list, part_coll_t * p); 
void dll_append(dll_t * list, part_coll_t * p);

dll_node_t * dll_at(dll_t * list, int index);

void dll_init_collisions(dll_t * list);

bool dll_extract_at(dll_t * list, int index, part_coll_t * p);
void dll_extract(dll_t * list, dll_node_t *node, part_coll_t * p);

void dll_delete_at(dll_t * list, int index);
void dll_delete(dll_t * list, dll_node_t *node); 

void dll_print(dll_t * list);
void dll_print_node(dll_node_t * list);

void dll_empty(dll_t * list);
dll_t * dll_copy_list(dll_t * rhs);

int dll_is_empty(dll_t * list);

#endif /* __DOUBLYLINKEDLIST_H */
