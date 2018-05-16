#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dll.h"

dll_t *dll_init(void)
{
    // Init our list
    dll_t * list = (dll_t *) malloc(sizeof *list);

    // Init the head
    list->head = (dll_node_t * ) malloc (sizeof *list->head);
    list->head->prev = NULL;

    // Init the tail
    list->tail = (dll_node_t * )  malloc (sizeof *list->tail);
    list->tail->next = NULL;

    // Head points to tail
    list->head->next = list->tail; 
    // Tail points back to head
    list->tail->prev = list->head; 

    // 0 elements at the beginning
	list->count = 0;

	return list;
}

void dll_destroy(dll_t * list)
{
    // Empty the list
	dll_empty(list);

    // Free special nodes and the actual list
    free(list->head);
    free(list->tail);
	free(list);
}

void dll_insert_after(dll_t * list, dll_node_t *node, pcord_t * p)
{
    dll_node_t * new_node = (dll_node_t *) malloc (sizeof *new_node);

    list->count++;

    new_node->prev = node;
    new_node->next = node->next;

    if (node->next == list->tail)
    {
        list->tail->prev = new_node;
    }
    else
    {
        node->next->prev = new_node;
    }

    node->next = new_node;
    new_node->p = p;
}

void dll_insert_before(dll_t * list, dll_node_t *node, pcord_t * p)
{
    dll_node_t * new_node = (dll_node_t *) malloc (sizeof *new_node);

    list->count++;
    new_node->prev = node->prev;
    new_node->next = node;

    if (node->prev == list->head)
    {
        list->head->next = new_node;
    }
    else
    {
        node->prev->next = new_node;
    }
    node->prev = new_node;
    new_node->p = p;
}

void dll_insert_beginning(dll_t * list, pcord_t * p)
{
    dll_insert_before(list, list->head->next, p);
}

void dll_prepend(dll_t * list, pcord_t * p)
{
    dll_insert_beginning(list, p);
}

void dll_insert_end(dll_t * list, pcord_t * p)
{
    dll_insert_after(list, list->tail->prev, p);
}

void dll_append(dll_t * list, pcord_t * p)
{
    dll_insert_end(list, p);
}

void dll_empty(dll_t * list)
{
    dll_node_t * current = list->head->next;
    while (current != list->tail)
	{
        current = dll_delete(list, current);
	}
    // Head points to tail
    list->head->next = list->tail; 
    // Tail points back to head
    list->tail->prev = list->head; 
    // And reset the count
	list->count = 0;
}

void dll_print(dll_t * list)
{
    printf("Size of list: %d\n", list->count);
    // Do stuff if something in the list
    if (!list->count) return;
    dll_node_t * node = list->head->next;
    int count = 0 ;
    while (node != list->tail)
    {
        pcord_t * ptr = node->p;
        printf("Particle %d: x=%.02f,y=%.02f,\tvx=%.02f,vy=%.02f.\n", count++,
                ptr->x, ptr->y, ptr->vx, ptr->vy);
        node = node->next;
    }
}

void dll_print_node(dll_node_t * node)
{
    if (!node)
    {
        printf("Not found.\n");
        return;
    }
    pcord_t * ptr = node->p;
    printf("Particle : x=%.02f,y=%.02f,\tvx=%.02f,vy=%.02f.\n",
            ptr->x, ptr->y, ptr->vx, ptr->vy);
}

void dll_append_list(dll_t * dst, dll_t * src)
{
    dll_node_t * current = src->head->next;
    while (current != src->tail)
    {
        // Allocate memory for the new particle
        pcord_t * p = (pcord_t *) malloc (sizeof *p);
        *p = *current->p;
        dll_append(dst, p);
        current = current->next;
    }
}

dll_t * dll_copy_list(dll_t * rhs)
{
    dll_t * list = dll_init();
    dll_node_t * nd = rhs->head->next;
    while (nd != rhs->tail)
    {
        dll_insert_end(list, nd->p);
        nd = nd->next;
    }
    return list;
}

int dll_is_empty(dll_t * list)
{
    return list->count == 0;
}

dll_node_t * dll_extract(dll_t * list, dll_node_t *node, pcord_t * p)
{
    p = node->p;
    return dll_delete(list, node);
}

dll_node_t * dll_extract_at(dll_t * list, int index, pcord_t * p)
{
    dll_node_t * node = dll_at(list, index);
    if (!node)
        return NULL;
    return dll_extract(list, node, p);
}

dll_node_t * dll_delete(dll_t * list, dll_node_t * node)
{
    dll_node_t * out;

    // Nothing to do here
    if (list->count == 0)
    {
        return NULL;
    }

    if (node->prev == list->head)
    {
        list->head->next = node->next;
    }
    else
    {
        node->prev->next = node->next;
    }
    if (node->next == list->tail)
    {
        list->tail->prev = node->prev;
    }
    else
    {
        node->next->prev = node->prev;
    }
    // Update the next node
    out = node->next;

    // Free memory
    free(node->p);
    free(node);

    // Decrease the count
    list->count--;

    return out;
}

dll_node_t * dll_delete_at(dll_t * list, int index)
{
    dll_node_t * node = dll_at(list, index);
    if (!node)
        return NULL;
    return dll_delete(list, node);
}

dll_node_t * dll_at(dll_t * list, int index)
{
    if (index < 0 || index >= list->count)
    {
        return NULL;
    }
    dll_node_t * node = list->head->next;
    int count = 0;
    while (node != list->tail)
    {
        if (count++ == index)
            return node;
        node = node->next;
    }
    return NULL;
}

pcord_t * dll_to_array(dll_t * list)
{
    pcord_t * outarray;
    if (list->count == 0)
       return NULL;
    outarray = (pcord_t * ) malloc (sizeof *outarray * list->count);
    int pos = 0;
    dll_node_t * node = list->head->next;
    while (node != list->tail)
    {
//         *(outarray + pos++) = *node->p;
        memcpy(outarray + pos++, node->p, sizeof(pcord_t));
        node = node->next;
    }
    return outarray;
}

dll_t * dll_from_array(pcord_t *  array, int count)
{
    dll_t * outlist = dll_init();
    for (int i = 0; i < count; ++i)
    {
        pcord_t * p = (pcord_t * ) malloc (sizeof *p);
        *p = *(array + i);
		dll_append(outlist, p);
	}
	return outlist;
}

bool dll_swap_nodes(dll_t * list, dll_node_t * node1, dll_node_t * node2)
{
    // Don't do this please...
    if (node1 == node2)
        return false;

    if (node1->next == node2)
    {
        node1->next = node2->next;
        node2->prev = node1->prev;

        if (node1->next != list->tail)
            node1->next->prev = node1;
        else
            list->tail->prev = node1;

        if (node2->prev != list->head)
            node2->prev->next = node2;
        else
            list->head->next = node2;

        node2->next = node1;
        node1->prev = node2;
    }
    else
    {
        dll_node_t * p = node2->prev;
        dll_node_t * n = node2->next;

        node2->prev = node1->prev;
        node2->next = node1->next;

        node1->prev = p;
        node1->next = n;

        if (node2->next != list->tail)
            node2->next->prev = node2;
        else
            list->tail->prev = node2;

        if (node2->prev != list->head)
            node2->prev->next = node2;
        else
            list->head->next = node2;

        if (node1->next != list->tail)
            node1->next->prev = node1;
        else
            list->tail->prev = node1;

        if (node1->prev != list->head)
            node1->prev->next = node1;
        else
            list->head->next = node1;

    }
    return true;
}

bool dll_swap(dll_t * list, int index1, int index2)
{
    if (index1 < 0 || index1 >= list->count || index2 < 0 || index2 >= list->count)
        return false;
    // Swap these nodes, give the indexes in order
    int min, max;
    if (index1 < index2)
    {
        min = index1;
        max = index2;
    }
    else
    {
        min = index2;
        max = index1;
    }
    dll_swap_nodes(list, dll_at(list, min), dll_at(list, max));
    return true;
}

