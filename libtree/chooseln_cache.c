#include<mathfunc.h>
#include <memalloc.h>
#include "chooseln_cache.h"

int chooseln_is_init2(struct chooseln_cache *cache)
{
	return cache->values ? 1 : 0;
}

int get_chooseln_cache_size2(struct chooseln_cache *cache)
{
	return cache->size;
}


void chooseln_cache_preset(struct chooseln_cache *cache, int maxFamilysize, int sfrom)
{
	int s, c, j;
	for (s = sfrom; s <= maxFamilysize; s++)
	{
		for (c = 0; c <= maxFamilysize; c++)
		{
			int m = MIN(s, c);
			double s_add_c = s + c;
			double s_add_c_sub_1 = s_add_c - 1;
			double s_sub_1 = s - 1;
			for (j = 0; j <= m; j++)
			{
				chooseln_get2(cache, s, j);
				chooseln_get2(cache, s_add_c_sub_1 - j, s_sub_1);
			}
		}
	}
}

void chooseln_cache_resize2(struct chooseln_cache *cache, int resize)
{
	if (cache->size >= resize) return;
	int i;
	if (cache->values)
	{
		cache->values = (double**)memory_realloc(cache->values, resize * 2, sizeof(double*));
	}
	else
	{
		cache->values = (double**)memory_new(resize * 2, sizeof(double*));
	}
	for (i = cache->size * 2; i < resize * 2; i++)
	{
		cache->values[i] = NULL;
	}
	int oldsize = cache->size;
	cache->size = resize;
	chooseln_cache_preset(cache, resize, oldsize + 1);
	fprintf(stderr, "** Cache resize: %d ==> %d\n", oldsize, resize);
}

void chooseln_cache_init2(struct chooseln_cache *cache, int size)
{
	cache->size = size;
	cache->values = (double**)memory_new(size * 2, sizeof(double*));
	int i;
	for (i = 0; i < size * 2; i++)
	{
		cache->values[i] = NULL;
	}
	chooseln_cache_preset(cache, size, 1);
}

void chooseln_cache_free2(struct chooseln_cache *cache)
{
	int i;
	for (i = 0; i < cache->size * 2; i++)
	{
		if (cache->values[i]) 
			memory_free(cache->values[i]);
		cache->values[i] = NULL;
	}
	memory_free(cache->values);
	cache->values = NULL;
}




