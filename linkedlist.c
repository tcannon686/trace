#include <stdlib.h>
#include <assert.h>

#include "linkedlist.h"



list_t *ListNew(void (*free)(void *), int (*cmp)(void *, void  *))
{
	list_t *ret = (list_t *)malloc(sizeof(list_t));
	ret->base_ptr = NULL;
	ret->top_ptr = NULL;
	ret->count = 0;
	ret->free = free;
	ret->cmp = cmp;
	return ret;
}

void ListFree(list_t *list)
{
	list_element_t *current = list->base_ptr;
	list_element_t *next_ptr;
	while(current != NULL)
	{
		next_ptr = current->next_ptr;
		if(list->free != NULL)
			list->free(current->data.data_ptr);
		free(current);
		current = next_ptr;
	}
	free(list);
}

void ListAppend(list_t *list_ptr, list_data_t data)
{
	list_element_t *element = (list_element_t *)malloc(sizeof(list_element_t));
	element->data = data;
	element->next_ptr = NULL;
	element->last_ptr = NULL;
	
	list_ptr->count ++;
	
	if(list_ptr->base_ptr == NULL)
	{
		list_ptr->base_ptr = element;
		list_ptr->top_ptr = element;
	}
	else
	{
		list_ptr->top_ptr->next_ptr = element;
		element->last_ptr = list_ptr->top_ptr;
		list_ptr->top_ptr = element;
	}
}

list_data_t ListGet(list_t *list_ptr, int index)
{
	assert(index < list_ptr->count);
	list_element_t *element = list_ptr->base_ptr;
	for(int i = 0; i < index; i ++)
	{
		element = element->next_ptr;
	}
	return element->data;
}

list_data_t ListGetLast(list_t *list_ptr)
{
	assert(list_ptr->count > 0);
	return list_ptr->top_ptr->data;
}

list_data_t ListGetFirst(list_t *list_ptr)
{
	assert(list_ptr->count > 0);
	return list_ptr->base_ptr->data;
}

void ListRemove(list_t *list_ptr, int index, int should_free)
{
	assert(index < list_ptr->count);
	list_ptr->count --;
	list_element_t *element = list_ptr->base_ptr;
	for(int i = 0; i < index; i ++)
	{
		element = element->next_ptr;
	}
	
	if(element->next_ptr != NULL)
		element->next_ptr->last_ptr = element->last_ptr;
	if(element->last_ptr != NULL)
		element->last_ptr->next_ptr = element->next_ptr;
	else
	{
		list_ptr->base_ptr = element->next_ptr;
	}
	
	if(should_free)
		list_ptr->free(element->data.data_ptr);
	
	free(element);
}

void ListRemoveFirst(list_t *list_ptr, int should_free)
{
	assert(list_ptr->count > 0);
	list_ptr->count --;
	
	list_ptr->base_ptr->next_ptr->last_ptr = NULL;
	if(should_free)
		list_ptr->free(list_ptr->base_ptr->data.data_ptr);
	
	list_element_t *old_ptr = list_ptr->base_ptr;
	list_ptr->base_ptr = list_ptr->base_ptr->next_ptr;
	free(old_ptr);
}


void ListRemoveLast(list_t *list_ptr, int should_free)
{
	assert(list_ptr->count > 0);
	list_ptr->count --;
	
	list_ptr->top_ptr->last_ptr->next_ptr = NULL;
	if(should_free)
		list_ptr->free(list_ptr->top_ptr->data.data_ptr);
	free(list_ptr->top_ptr);
}

int ListIndexOf(list_t *list_ptr, list_data_t data)
{
	assert(list_ptr->cmp != NULL);
	
	list_element_t *element = list_ptr->base_ptr;
	int i = 0;
	while(element != NULL)
	{
		if(list_ptr->cmp(element->data.data_ptr, data.data_ptr) == 0)
			return i;
		i ++;
	}
	return -1;
}

