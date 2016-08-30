#ifndef __GMATRIX_H__
#define __GMATRIX_H__

typedef struct 
{
	int dim;		
	int* size;
	int* cumsize;
	int datasize;
	int num_elements;
	void** data;
}GMatrix;

typedef GMatrix* pGMatrix;

extern pGMatrix gmatrix_new(int dim, int* size, int datasize );
extern void gmatrix_free(pGMatrix pgm);
extern int gmatrix_index(pGMatrix pgm, ... );
extern void* gmatrix_get(pGMatrix pgm, ... );
extern void gmatrix_set(pGMatrix pgm, void* psrc, ... );
extern pGMatrix gmatrix_double_new(int dim, int* size);
extern double gmatrix_double_get(pGMatrix pgm, ... );
extern void gmatrix_double_set(pGMatrix pgm, double d, ... );
extern void gmatrix_double_set_with_index(pGMatrix pgm, double d, int idx );
extern void gmatrix_dim_index(pGMatrix pgm, int idx, int* didx );
extern double gmatrix_double_get_with_index(pGMatrix pgm, int idx );

#endif
