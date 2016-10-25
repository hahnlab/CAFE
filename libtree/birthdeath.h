#ifndef __BIRTHDEATH_H__
#define __BIRTHDEATH_H__

#include "hashtable.h"
#include "chooseln_cache.h"

typedef struct
{
	int branchlength;
	double lambda;
	double mu;
	int maxFamilysize;
}BirthDeathCacheAttrib;

typedef BirthDeathCacheAttrib* pBirthDeathCacheAttrib;

typedef struct
{
	//pArrayList list;
    hash_table_t* table;
	int maxFamilysize;
	//int base_bl;
}BirthDeathCacheArray;
typedef BirthDeathCacheArray* pBirthDeathCacheArray;

typedef struct 
{
	BirthDeathCacheAttrib attrib;
	double** matrix;	
}BirthDeathCache;
typedef BirthDeathCache* pBirthDeathCache;

extern void birthdeath_cache_array_free(pBirthDeathCacheArray pbdc_array);
extern double birthdeath_likelihood_with_s_c(int s, int c, double branchlength, double lambda, double mu, struct chooseln_cache *cache);
extern pBirthDeathCache eq_birthdeath_cache_new( double branchlength, double lambda, int maxFamilysize );
extern pBirthDeathCache birthdeath_cache_new( double branchlength, double lambda, double mu, int maxFamilysize );
double** eq_birthdeath_cache_get_matrix(pBirthDeathCacheArray pbdc_array, double branchlength, double lambda );
double** birthdeath_cache_get_matrix(pBirthDeathCacheArray pbdc_array, double branchlength, double lambda, double mu );
extern void thread_run(int numthreads, void* (*run)(void*), void* param, int size );
extern pBirthDeathCacheArray birthdeath_cache_array_new_with_list_thread(int* bl, int size, int maxFamilysize, double lambda, int numthreads );
double birthdeath_rate_with_log_alpha(int s, int c, double log_alpha, double coeff, struct chooseln_cache *cache);

/**
* \brief A cache of values of chooseln
*
* Chooseln evaluates the natural logarithm of Gamma(n+1)/(Gamma(k+1)*Gamma(n-k+1))
* The cache holds values for integer values of n and k. It does not appear to be 
* threadsafe.
*/
extern int chooseln_is_init();
extern int get_chooseln_cache_size();
extern void chooseln_cache_init(int size);
extern void chooseln_cache_resize(int resize);
extern void chooseln_cache_free();
#endif
