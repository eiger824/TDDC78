#ifndef __DOUBLYLINKEDLIST_H
#define __DOUBLYLINKEDLIST_H

#include "definitions.h"

typedef struct dll_node_type
{
    pcord_t * p;
    struct dll_node_type *next;
    struct dll_node_type *prev;

} dll_node_t;

typedef struct dll_type
{
    dll_node_t * head;
    dll_node_t * tail;

    int count;

} dll_t;

dll_t *dll_init(void);
void dll_destroy(dll_t * list);

void dll_insert_after(dll_t * list, dll_node_t *node, pcord_t * p);
void dll_insert_before(dll_t * list, dll_node_t *node, pcord_t * p);

void dll_insert_beginning(dll_t * list, pcord_t * p); 
void dll_prepend(dll_t * list, pcord_t * p);

void dll_insert_end(dll_t * list, pcord_t * p); 
void dll_append(dll_t * list, pcord_t * p);

dll_node_t * dll_at(dll_t * list, int index);

dll_node_t * dll_extract_at(dll_t * list, int index, pcord_t * p);
dll_node_t * dll_extract(dll_t * list, dll_node_t *node, pcord_t * p);

dll_node_t * dll_delete_at(dll_t * list, int index);
dll_node_t * dll_delete(dll_t * list, dll_node_t * node);

void dll_print(dll_t * list);
void dll_print_node(dll_node_t * node);

void dll_empty(dll_t * list);
dll_t * dll_copy_list(dll_t * rhs);
void dll_append_list(dll_t * dst, dll_t * src);

int dll_is_empty(dll_t * list);

pcord_t * dll_to_array(dll_t * list);
dll_t * dll_from_array(pcord_t * array, int count);

bool dll_swap_nodes(dll_t * list, dll_node_t * node1, dll_node_t * node2);
bool dll_swap(dll_t * list, int index1, int index2);

#endif /* __DOUBLYLINKEDLIST_H */

