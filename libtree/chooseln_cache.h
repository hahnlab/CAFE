#ifndef CHOOSELN_CACHE_H_2C3C614A_C14D_477A_B8A0_30AE50238D9E
#define CHOOSELN_CACHE_H_2C3C614A_C14D_477A_B8A0_30AE50238D9E
/**
* \brief A cache of values of chooseln
*
* Chooseln evaluates the natural logarithm of Gamma(n+1)/(Gamma(k+1)*Gamma(n-k+1))
* The cache holds values for integer values of n and k. It does not appear to be
* threadsafe.
*/

#include <mathfunc.h>
#include <memalloc.h>

struct chooseln_cache {
	double** values;
	int size;
};

int chooseln_is_init2(struct chooseln_cache *cache);
int get_chooseln_cache_size2(struct chooseln_cache *cache);
void chooseln_cache_init2(struct chooseln_cache *cache, int size);
void chooseln_cache_resize2(struct chooseln_cache *cache, int resize);
void chooseln_cache_free2(struct chooseln_cache *cache);

static inline double chooseln_get2(struct chooseln_cache *cache, int n, int x)
{
	if (cache->values[n] && cache->values[n][x] >= 0) return cache->values[n][x];
	if (cache->values[n] == NULL)
	{
		int i;
		cache->values[n] = (double*)memory_new(cache->size + 1, sizeof(double));
		for (i = 0; i <= cache->size; i++) 
			cache->values[n][i] = -1.0;
	}
	cache->values[n][x] = chooseln(n, x);
	return cache->values[n][x];
}

#endif
