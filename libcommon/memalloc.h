#ifndef __MEMALLOC_H__
#define __MEMALLOC_H__

#include<string.h>
#include<errno.h>
#include<stdarg.h>
#include<stdlib.h>
#include<stdio.h>


typedef void (*func_memory_free)(void* data);

extern void* memory_new(size_t count, size_t size );
extern void* memory_new_with_init(int count, int size, void* init);
extern void* memory_realloc_with_data(void** data, int old, int newsize, int size, func_memory_free func, int elementsize );
extern void* memory_realloc(void* data, int count, int size );
extern void  memory_free(void* data);
extern void** memory_new_2dim(int row, int col, int size );
extern void memory_free_2dim(void** data, int row, int col, func_memory_free func );
extern void** memory_copy_2dim(void** dst, void** src, int row, int col, int size );

#endif 
