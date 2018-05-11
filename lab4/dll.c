#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dll.h"

dll_t *dll_init(void)
{
    dll_t * list = (dll_t *) malloc(sizeof(dll_t));
    list->first = NULL;
    list->last = NULL;
    list->count = 0;
    return list;
}

void dll_destroy(dll_t * list)
{
    dll_empty(list);
    free(list);
}

void dll_insert_after(dll_t * list, dll_node_t *node, pcord_t  p)
{
    dll_node_t *new_node = (dll_node_t *) malloc(sizeof(dll_node_t));
    list->count++;
    new_node->prev = node;
    new_node->next = node->next;
    if (node->next == NULL)
    {
        list->last = new_node;
    }
    else
    {
        node->next->prev = new_node;
    }
    node->next = new_node;
    new_node->p = p;
}

void dll_insert_before(dll_t * list, dll_node_t *node, pcord_t  p)
{
    dll_node_t *new_node = (dll_node_t *) malloc(sizeof(dll_node_t));
    list->count++;
    new_node->prev = node->prev;
    new_node->next = node;
    if (node->prev == NULL)
    {
        list->first = new_node;
    }
    else
    {
        node->prev->next = new_node;
    }
    node->prev = new_node;
    new_node->p = p;
}

void dll_insert_beginning(dll_t * list, pcord_t  p)
{
    if (list->first == NULL)
    {
        dll_node_t *new_node = (dll_node_t *) malloc(sizeof(dll_node_t));
        list->count++;
        list->first = new_node;
        list->last = new_node;
        new_node->next = NULL;
        new_node->prev = NULL;
        new_node->p = p;
    }
    else
    {
        dll_insert_before(list, list->first, p);
    }
}

void dll_prepend(dll_t * list, pcord_t  p)
{
    dll_insert_beginning(list, p);
}

void dll_insert_end(dll_t * list, pcord_t  p)
{
    if (list->last == NULL)
    {
        dll_insert_beginning(list, p);
    }
    else
    {
        dll_insert_after(list, list->last, p);
    }
}

void dll_append(dll_t * list, pcord_t  p)
{
    dll_insert_end(list, p);
}

void dll_empty(dll_t * list)
{
    dll_node_t *temp_first_node;
    while (list->first != NULL)
    {
        temp_first_node = list->first;
        list->first = list->first->next;
        free(temp_first_node);
    }
    list->count = 0;
    list->last = NULL;
}

void dll_print(dll_t * list)
{
    printf("Size of list: %d\n", list->count);
    // Do stuff if something in the list
    if (!list->count) return;
    dll_node_t * node = list->first;
    int count = 0 ;
    while (node != list->last->next)
    {
        pcord_t ptr = node->p;
        printf("Particle %d: x=%.02f,y=%.02f,\tvx=%.02f,vy=%.02f.\n", count++,
                ptr.x, ptr.y, ptr.vx, ptr.vy);
        node = node->next;
    }
}

void dll_print_node(dll_node_t * list)
{
    if (!list)
    {
        printf("Not found.\n");
        return;
    }
    pcord_t ptr = list->p;
    printf("Particle : x=%.02f,y=%.02f,\tvx=%.02f,vy=%.02f.\n",
            ptr.x, ptr.y, ptr.vx, ptr.vy);
}

dll_t * dll_copy_list(dll_t * rhs)
{
    dll_t * list = dll_init();
    dll_node_t * nd = rhs->first;
    while (nd != rhs->last->next)
    {
        dll_insert_end(list, nd->p);
        nd = nd->next;
    }
    return list;
}

int dll_is_empty(dll_t * list)
{
    return !list->count;
}

bool dll_extract(dll_t * list, dll_node_t *node, pcord_t * p)
{
    *p = node->p;
    return dll_delete(list, node);
}

bool dll_extract_at(dll_t * list, int index, pcord_t * p)
{
    dll_node_t * node = dll_at(list, index);
    if (!node)
        return false;
    return dll_extract(list, node, p);
}

/* Return true if there are still elements on the list */
bool dll_delete(dll_t * list, dll_node_t *node)
{
    if (node->prev == NULL)
    {
        list->first = node->next;
    }
    else
    {
        node->prev->next = node->next;
    }
    if (node->next == NULL)
    {
        list->last = node->prev;
    }
    else
    {
        node->next->prev = node->prev;
    }
    free(node);
    return (!--list->count ? false: true);
}

bool dll_delete_at(dll_t * list, int index)
{
    dll_node_t * node = dll_at(list, index);
    if (!node)
        return false;
    return dll_delete(list, node);
}

dll_node_t * dll_at(dll_t * list, int index)
{
    if (index < 0 || index >= list->count)
    {
        return NULL;
    }
    dll_node_t * node = list->first;
    int count = 0;
    while (node != list->last->next)
    {
        if (count++ == index)
            return node;
        node = node->next;
    }
    return NULL;
}

pcord_t * dll_to_array(dll_t * list)
{
    pcord_t * outarray = (pcord_t * ) malloc (sizeof *outarray * list->count);
    int pos = 0;
    dll_node_t * node = list->first;
    while (node != list->last->next)
    {
        *(outarray + pos++) = node->p;
        node = node->next;
    }
    return outarray;
}

dll_t * dll_from_array(pcord_t *  array, int count)
{
    dll_t * outlist = dll_init();
    for (int i = 0; i < count; ++i)
    {
        pcord_t p = *(array + i);
		dll_append(outlist, p);
	}
	return outlist;
}

