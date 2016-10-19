#ifndef __UTILS_H__
#define __UTILS_H__

#include<memalloc.h>
#include<stdint.h>

typedef void (*freefunc)(void*);
/********************************************************************
 * ArrayList
 ********************************************************************/
typedef struct
{
	void** array;
	int	size;
	int remain;
	int step;
}ArrayList;
typedef ArrayList* pArrayList;

extern pArrayList arraylist_new(int step);
extern pArrayList arraylist_add(pArrayList pal, void* data);
extern void* arraylist_get(pArrayList pal, int idx);
extern void arraylist_free(pArrayList pal, freefunc datafree );
extern void arraylist_trim(pArrayList pal);
extern void arraylist_clear(pArrayList pal);
extern void arraylist_sort(pArrayList pal, int (*compar)(const void*, const void*) );
extern void arraylist_shuffle(pArrayList pal);


/********************************************************************
 * LinkedLIstItem 
 ********************************************************************/
typedef struct tagLinkedListItem* pLinkedListItem;
typedef struct tagLinkedListItem
{
	void* 			data;
	pLinkedListItem	prev;
	pLinkedListItem	next;
}LinkedListItem;

extern pLinkedListItem linkedlistitem_new(void* data);
extern void linkedlistitem_free(pLinkedListItem pitem);

/********************************************************************
 * Vector
 ********************************************************************/
typedef struct tagVector* pVector;
typedef struct tagVector
{
	int 			size;
	pLinkedListItem	head;
	pLinkedListItem	tail;
	pLinkedListItem cur;
}Vector;

extern pVector vector_new();
extern void vector_free(pVector pvec, freefunc func);
extern void vector_add(pVector pvec, void* data);
extern void* vector_get(pVector pvec, int idx);
extern void vector_rewind(pVector pvec);
extern void* vector_next(pVector pvec);
extern void vector_sort(pVector pvec, int (*cmp)(const void*,const void*) );
extern void* vector_get_by_cmp( pVector pvec, void* data, size_t (*cmp)(const void*,const void*) ) ;
extern void vector_remove_by_data(pVector pvec, void* data);
extern void vector_dereference_by_data(pVector pvec, void* data);
extern pArrayList vector_to_arraylist(pVector pvec);

typedef Vector Stack;
typedef pVector pStack;

extern pStack stack_new();
extern void stack_free(pStack pstack);
extern void stack_push(pStack pstack, void* data);
extern void* stack_pop(pStack pstack);
extern int stack_is_empty(pStack pstack);


/********************************************************************
 * Etc
 ********************************************************************/
extern void print_error(char* file, char* function, int line, char* message, ... );
extern int __cmp_int(const void* a, const void* b);
extern int __cmp_double(const void* a, const void* b);

#endif
