#include<stdlib.h>
#include<stdio.h>
#include<stdarg.h>
#include<errno.h>
#include<time.h>
#include "utils.h"
#include "memalloc.h"

int __cmp_int(const void* a, const void* b)
{
	return *(int*)a - *(int*)b;
}

int __cmp_double(const void* a, const void* b) 
{           
	if ( *(double*)a > *(double*)b ) return 1;
	if ( *(double*)a == *(double*)b ) return 0;
	return -1;
}   

void print_error(char* file, char* function, int line, char* message, ... )
{
	va_list ap;
	va_start(ap, message);
	char buf[BUFSIZ];
	vsprintf(buf,message,ap);
	va_end(ap);
	fprintf(stderr,"%s:%s:%d: %s\n", file,function,line, buf );
	exit(errno);
}

/********************************************************************
 * ArrayList
 ********************************************************************/

#define ARRAYLIST_DEFAULT_STEP 1024

pArrayList arraylist_new(int step)
{
	if ( step <= 0 )
	{
		step = ARRAYLIST_DEFAULT_STEP;
	}
	pArrayList pal = (pArrayList)memory_new(1, sizeof(ArrayList));
	pal->size = 0;
	pal->remain = step;
	pal->array = (void**)memory_new(step, sizeof(void*));
	pal->step = step;
	return pal;
}

pArrayList arraylist_add(pArrayList pal, void* data)
{
	pal->remain--;
	if ( pal->remain == 0 )
	{
		pal->array = (void**)realloc(pal->array, sizeof(void*)*(pal->size + pal->step) );
		if ( pal->array == NULL )
		{
			print_error(__FILE__,(char*)__FUNCTION__,__LINE__,
					    "memory reallocation of ArrayList->array");		
		}
		pal->remain = pal->step;
	}
	pal->array[pal->size++] = data;
	return pal;
}

void arraylist_clear(pArrayList pal)
{
	pal->remain += pal->size;
	pal->size = 0;
}

void arraylist_trim(pArrayList pal)
{
	pal->array = (void**)realloc(pal->array,sizeof(void*)*pal->size);
	pal->remain = 0;
}

void arraylist_sort(pArrayList pal, int (*compar)(const void*, const void*) )
{
	qsort(pal->array, pal->size, sizeof(void*), compar );
}

void* arraylist_get(pArrayList pal, int idx)
{
	if ( idx >= pal->size ) return NULL;	
	return 	pal->array[idx];
}

void arraylist_free(pArrayList pal, freefunc datafree )
{
	if ( datafree )
	{
		int i = 0;
		for ( i = 0 ; i < pal->size ; i++ )
		{
			if ( pal->array[i] ) datafree(pal->array[i]);
            pal->array[i] = NULL;
		}
	}
	memory_free(pal->array);
	pal->array = NULL;
	memory_free(pal);
	pal = NULL;
}


void arraylist_shuffle(pArrayList pal) 
{
	int i;
	for(i = 0; i < pal->size-1; i++) {
		int c = rand() / (RAND_MAX/(pal->size-i) + 1); 
		void* t = pal->array[i]; pal->array[i] = pal->array[i+c]; pal->array[i+c] = t;	/* swap */
	}
}


/********************************************************************
 * LinkedLIstItem 
 ********************************************************************/
pLinkedListItem linkedlistitem_new(void* data)
{
	pLinkedListItem pitem = (pLinkedListItem)memory_new(1,sizeof(LinkedListItem));
	pitem->data = data;
	pitem->next = NULL;
	pitem->prev = NULL;
	return pitem;
}

void linkedlistitem_free(pLinkedListItem pitem)
{
	if ( pitem->prev )
	{
		pitem->prev->next = pitem->next;
	}
	if ( pitem->next )
	{
		pitem->next->prev = pitem->prev;
	}
	memory_free(pitem);
	pitem = NULL;
}

void linkedlistitem_dereference(pLinkedListItem pitem)
{
	if ( pitem->prev )
	{
		pitem->prev->next = pitem->next;
	}
	if ( pitem->next )
	{
		pitem->next->prev = pitem->prev;
	}
}

/********************************************************************
 * Vector
 ********************************************************************/
pVector vector_new()
{
	pVector pvec = (pVector)memory_new(1,sizeof(Vector));
	pvec->size = 0;
	pvec->head = NULL;
	pvec->tail = NULL;
	return pvec;
}

void vector_free(pVector pvec, freefunc datafree)
{
	int i = 0;
	pLinkedListItem cur = pvec->head;
	pLinkedListItem tmp;
	for( i = 0 ; i < pvec->size; i++ )
	{
		tmp = cur->next;
		if ( datafree ) datafree(cur->data);
		memory_free(cur);
		cur = tmp;
	}
	memory_free(pvec);
	pvec = NULL;
}

void vector_add(pVector pvec, void* data)
{
	pvec->size++;
	pLinkedListItem pitem = linkedlistitem_new(data);
	if ( pvec->head == NULL )
	{
		pvec->head = pvec->tail = pitem;
	}
	else
	{
		pitem->prev = pvec->tail;
		pvec->tail->next = pitem;
		pvec->tail = pitem;
	}
}

void* vector_get(pVector pvec, int idx)
{
	int i = 0;
	pLinkedListItem cur;
	if ( idx >= pvec->size ) return NULL;
	cur = pvec->head;
	for ( i = 0 ; i < idx ; i++ )
	{
		cur = cur->next;
	}
	return cur->data;
}

void* vector_get_by_cmp( pVector pvec, void* data, size_t (*cmp)(const void*,const void*) ) 
{
	int i = 0;
	pLinkedListItem cur;
	cur = pvec->head;
	for ( i = 0 ; i < pvec->size ; i++ )
	{
		if ( cmp( cur->data, data) == 0 )
		{
			return cur->data;
		}
		cur = cur->next;
	}
	return NULL;
}

void vector_sort(pVector pvec, int (*cmp)(const void*,const void*) )
{
	pArrayList pal = vector_to_arraylist(pvec);
	arraylist_sort(pal,cmp);
	int i;
	pLinkedListItem cur = pvec->head;
	for ( i = 0 ; i < pal->size ; i++ )
	{
		cur->data = pal->array[i];
		cur = cur->next;
	}
	arraylist_free(pal,NULL);
}

void vector_remove_by_data(pVector pvec, void* data)
{
	int i;
	pLinkedListItem cur;
	if ( pvec->size == 0 ) return;
	cur = pvec->head;
	for ( i = 0 ; i < pvec->size ; i++ )
	{
		if ( cur->data == data ) break;
	}
	pvec->size--;
	if (data == pvec->head->data) {
		pvec->head = pvec->head->next;		
	}
	if (data == pvec->tail->data) {
		pvec->tail = pvec->tail->prev;
	}
	linkedlistitem_free(cur);
}

void vector_dereference_by_data(pVector pvec, void* data)
{
	int i;
	pLinkedListItem cur;
	if ( pvec->size == 0 ) return;
	cur = pvec->head;
	for ( i = 0 ; i < pvec->size ; i++ )
	{
		if ( cur->data == data ) break;
	}
	pvec->size--;
	if (data == pvec->head->data) {
		pvec->head = pvec->head->next;		
	}
	if (data == pvec->tail->data) {
		pvec->tail = pvec->tail->prev;
	}
	linkedlistitem_dereference(cur);
}

pArrayList vector_to_arraylist(pVector pvec)
{
	if ( pvec->size == 0 ) return NULL;
	pArrayList pal = arraylist_new(pvec->size);
	pLinkedListItem cur;
	cur = pvec->head;
	int i;
	for ( i = 0 ; i < pvec->size ; i++ )
	{
		arraylist_add(pal,cur->data);
		cur = cur->next;
	}
	return pal;
}

void vector_rewind(pVector pvec)
{
	pvec->cur = pvec->head;		
}

void* vector_next(pVector pvec)
{
	if ( pvec->cur == NULL ) return NULL;
	void* data = pvec->cur->data;
	pvec->cur = pvec->cur->next;
	return data;
}

/********************************************************************
 * Stack
 ********************************************************************/
pStack stack_new()
{
	return vector_new();
}

void stack_free(pStack pstack)
{
	vector_free(pstack, NULL);
}

void stack_push(pStack pstack, void* data)
{
	pstack->size++;
	pLinkedListItem pitem = linkedlistitem_new(data);
	if ( pstack->head )
	{
		pitem->next = pstack->head;	
		pstack->head->prev = pitem;
		pstack->head = pitem;
	}
	else
	{
		pstack->head = pstack->tail = pitem;
	}
}

void* stack_pop(pStack pstack)
{
	if ( pstack->size == 0 )
	{
		return NULL;
	}
	pstack->size--;
	pLinkedListItem head = pstack->head;
	pstack->head = head->next;
	void* data  = head->data;
	linkedlistitem_free(head);
	if ( pstack->size == 0 )
	{
		pstack->head = pstack->tail = NULL;
	}
	return data;	
}

int stack_has_items(pStack pstack)
{
	return pstack->size;
}

