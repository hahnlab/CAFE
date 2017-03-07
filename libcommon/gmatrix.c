#include<stdarg.h>

#include "gmatrix.h"
#include "utils.h"
#include "memalloc.h"

pGMatrix gmatrix_new(int dim, int* size, int datasize )
{
	int i;
	pGMatrix pgm = (pGMatrix)memory_new(1,sizeof(GMatrix));
	pgm->dim = dim;
	pgm->size = (int*)memory_new(dim,sizeof(int));
	pgm->cumsize = (int*)memory_new(dim,sizeof(int));
	memcpy(pgm->size,size,dim*sizeof(int));
	pgm->cumsize[dim-1] = 1;
	int total = 1;
	for ( i = 0 ; i < dim ; i++ ) total *= size[i];
	for ( i = dim - 1 ; i > 0 ; i-- ) pgm->cumsize[i-1] = pgm->cumsize[i] * size[i];
	pgm->data = (void**)memory_new( total, datasize );
	pgm->num_elements = total;
	return pgm;
}

void gmatrix_free(pGMatrix pgm)
{
	memory_free(pgm->data);
	pgm->data = NULL;
	memory_free(pgm->size);
	pgm->size = NULL;
	memory_free(pgm->cumsize);
	pgm->cumsize = NULL;
	memory_free(pgm);	
	pgm = NULL;
}

int gmatrix_vindex(pGMatrix pgm, va_list ap1)
{
	int i, idx;
  va_list ap;
  va_copy(ap, ap1);
	idx = 0;
	for ( i = 0; i < pgm->dim - 1 ; i++ )
	{
		idx += pgm->cumsize[i] * va_arg(ap,int);			
	}
	idx += va_arg(ap,int);
	va_end(ap);
	return idx;
}

void gmatrix_dim_index(pGMatrix pgm, int idx, int* didx )
{
	int i;
	for ( i = 0 ; i < pgm->dim ; i++ )
	{
		didx[i] = idx / pgm->cumsize[i];	
		idx -= didx[i] * pgm->cumsize[i];
	}
}

int gmatrix_index(pGMatrix pgm, ... )
{
	va_list ap;	
	va_start(ap,pgm);
	int idx = gmatrix_vindex(pgm,ap);	
	va_end(ap);
	return idx;
}

void* gmatrix_vget(pGMatrix pgm, va_list ap1)
{
  va_list ap;
  va_copy(ap, ap1);
	uintptr_t idx = gmatrix_vindex(pgm,ap);	
  va_end(ap);
	return (void*)(((uintptr_t)pgm->data)+idx*pgm->datasize);
}

void* gmatrix_get(pGMatrix pgm, ... )
{
	va_list ap;	
	va_start(ap,pgm);
	void* data = gmatrix_vget(pgm,ap);
	va_end(ap);
	return data;
}

void gmatrix_set(pGMatrix pgm, void* psrc, ... )
{
	va_list ap;	
	va_start(ap,psrc);
	void* data = gmatrix_vget(pgm, ap);
	va_end(ap);
	memcpy( data, psrc, pgm->datasize );
}

pGMatrix gmatrix_double_new(int dim, int* size)
{
	return gmatrix_new(dim,size, sizeof(double));
}

double gmatrix_double_get(pGMatrix pgm, ... )
{
	va_list ap;	
	va_start(ap,pgm);
	void* data = gmatrix_vget(pgm,ap);
	va_end(ap);
	return *(double*)data;
}

void gmatrix_double_set(pGMatrix pgm, double d, ... )
{
	va_list ap;	
	va_start(ap,d);
	void* data = gmatrix_vget(pgm, ap);
	va_end(ap);
	*(double*)data = d;
}

void gmatrix_double_set_with_index(pGMatrix pgm, double d, int idx )
{
	((double*)pgm->data)[idx] = d;
}

double gmatrix_double_get_with_index(pGMatrix pgm, int idx )
{
	return (double)((double*)pgm->data)[idx];
}
