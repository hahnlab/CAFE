#include "../config.h"
#include <assert.h>

#include "birthdeath.h"
#include<utils.h>
#include<stdlib.h>
#include<stdio.h>
#include<family.h>
#include<float.h>
#include<mathfunc.h>
#include "chooseln_cache.h"

#ifdef HAVE_BLAS
#ifdef HAVE_ATLAS
#include "cblas.h"
#else
#include "mkl.h"
#endif
#endif

/*
P(X(t) = c | X(0) = s)  = \sum_{j=0}^{\min(s,c)} \binom{s}{j}\binom{s+c-j-1}{s-1}
                                                 \alpha^{s+c-2j}(1-2\alpha)^j
 */

struct chooseln_cache cache = { 0,0 };

int chooseln_is_init()
{
	return chooseln_is_init2(&cache);
}

int get_chooseln_cache_size() 
{ 
	return get_chooseln_cache_size2(&cache);
}

double chooseln_get(int n, int x)
{
	return chooseln_get2(&cache, n, x);
}

void chooseln_cache_resize(int resize)
{
	chooseln_cache_resize2(&cache, resize);
}

void chooseln_cache_init(int size)
{
	chooseln_cache_init2(&cache, size);
}

void chooseln_cache_free()
{
	chooseln_cache_free2(&cache);
}




double birthdeath_rate_with_log_alpha_beta(int s, int c, double log_alpha, double log_beta, double log_coeff, struct chooseln_cache *cache)
{
	assert(cache->values != 0);
	int j;
	int m = MIN(c,s);
	double t, p = 0;
	int s_add_c = s + c;
	int s_add_c_sub_1 = s_add_c - 1;
	int s_sub_1 = s - 1;
	for ( j = 0, p = 0 ; j <= m ; j++ )
	{
		t = cache->values[s][j] + cache->values[s_add_c_sub_1-j][s_sub_1] + (s-j)*log_alpha + (c-j)*log_beta + j*log_coeff;
		//t = chooseln_get(s, j) + chooseln_get(s_add_c_sub_1-j,s_sub_1) + (s-j)*log_alpha + (c-j)*log_beta + j*log_coeff;
		p += exp(t);
	}
	return MAX(MIN(p,1),0);
}

double birthdeath_rate_with_log_alpha(int s, int c, double log_alpha, double coeff, struct chooseln_cache *cc )
{
	if (cc == NULL)
		cc = &cache;

	assert(cc->values != 0);
	int m = MIN(c,s);

	double lastterm = 1;
	double p = 0.0;
	int s_add_c = s + c;
	int s_add_c_sub_1 = s_add_c - 1;
	int s_sub_1 = s - 1;
	for (int j = 0; j <= m ; j++ )
	{
		double t = chooseln_get2(cc, s, j) + chooseln_get2(cc, s_add_c_sub_1-j,s_sub_1) + (s_add_c-2*j)*log_alpha;
		p += (exp(t) * lastterm);
		lastterm *= coeff;
	}

	return MAX(MIN(p,1),0);
}

/**
* \brief Calculates the probability of transitioning from root_family_size to family_size
*
* Given the branch length and the expected change rate lammbda
*/double birthdeath_likelihood_with_s_c(int root_family_size, int family_size, double branchlength, double lambda, double mu, struct chooseln_cache *cache)
{	
	double alpha, coeff, beta=0;
	double denominator, numerator = 0;	
	
	if (root_family_size == 0) {
		if (family_size == 0) {
			return 1;
		}
		else {
			return 0;
		}
	}
	else {
		if ((mu < 0) || (lambda == mu)) {	// this is worrisome. Comparing two doubles is not an exact process
			alpha = lambda*branchlength/(1+lambda*branchlength);
			coeff = 1 - 2*alpha;
			if (coeff <= 0) {
				printf("Negative coefficient found (lambda %f branchlength %f)\n", lambda, branchlength);
				return 0;
			}
			else {
				double result = birthdeath_rate_with_log_alpha(root_family_size, family_size,log(alpha),coeff, cache);
				return result;
			}
		}
		else {
			denominator = lambda*(exp((lambda-mu)*branchlength))-mu;
			numerator = exp((lambda-mu)*branchlength)-1;
			alpha = (mu*numerator)/denominator;
			beta = (lambda*numerator)/denominator;
			coeff = 1 - alpha - beta;
			if (coeff <= 0) {
				return 0;
			}
			else {
				return birthdeath_rate_with_log_alpha_beta(root_family_size, family_size, log(alpha),log(beta),log(coeff), cache);
			}
		}
	}
}

void square_matrix_init(struct square_matrix* matrix, int sz)
{
	matrix->values = (double*)memory_new(sz * sz, sizeof(double));
	matrix->size = sz;
}

void square_matrix_set(struct square_matrix* matrix, int x, int y, double val)
{
	assert(x < matrix->size);
	assert(y < matrix->size);
	matrix->values[x*matrix->size+y] = val;
}

void square_matrix_delete(struct square_matrix* matrix)
{
	memory_free((void*)matrix->values);
}

void square_matrix_resize(struct square_matrix* matrix, int new_size)
{
	int n = new_size < matrix->size ? new_size : matrix->size;
	double *new_values = memory_new(new_size*new_size, sizeof(double*));
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			new_values[i*new_size + j] = matrix->values[i*matrix->size + j];
	memory_free(matrix->values);
	matrix->values = new_values;
	matrix->size = new_size;
}

void square_matrix_print(struct square_matrix* matrix)
{
  for (int s = 0; s < matrix->size; s++)
  {
    for (int c = 0; c < matrix->size; c++)
    {
      printf("%e ", square_matrix_get(matrix, s, c));
    }
    printf("\n");
  }
}

void square_matrix_multiply(struct square_matrix* matrix, double *vector, int row_start, int row_end, int col_start, int col_end, double *result)
{
#ifdef HAVE_BLAS
  double alpha = 1.0, beta = 0.;
  int m = row_end - row_start + 1;
  int k = col_end - col_start + 1;
  int n = 1;
  double *sub = matrix->values + row_start*matrix->size + col_start; 
  cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,m,n,k,alpha,sub,matrix->size,vector,n,beta,result,n);
#else
  for (int s = row_start, i = 0; s <= row_end; s++, i++)
  {
    result[i] = 0;
    for (int c = col_start, j = 0; c <= col_end; c++, j++)
    {
      result[i] += square_matrix_get(matrix, s, c) * vector[j];
    }
  }
#endif
}

void init_zero_matrix(struct square_matrix *matrix)
{
	for (int s = 1; s < matrix->size; s++)
	{
		for (int c = 0; c < matrix->size; c++)
		{
			square_matrix_set(matrix, s, c, 0);
		}
	}
}

void init_identity_matrix(struct square_matrix* matrix)
{
	for (int s = 1; s < matrix->size; s++)
	{
		for (int c = 0; c < matrix->size; c++)
		{
			if (s == c) {
				square_matrix_set(matrix, s, c, 1);
			}
			else {
				square_matrix_set(matrix, s, c, 0);
			}
		}
	}
}

void init_matrix(struct square_matrix* matrix, double coeff)
{
	for (int c = 1; c < matrix->size; c++)
	{
		square_matrix_set(matrix, 0, c, 0);
	}
	if (coeff <= 0)
	{
		init_zero_matrix(matrix);
	}
	else if (coeff == 1)
	{
		init_identity_matrix(matrix);
	}
}

// THE FUNCTION!!!!!!!!!!
// must add mu to calculate the transition probability 
/**
	returns a structure representing a matrix of precalculated transition probabilites from one family
	size to another, given the specified values of branch length, lambda, and mu. Values are calculated
	from 0 to the given maxFamilySize. All relevant values are stored in the structure.

	lambda is the probability of both gene gain and loss per gene per unit time in the phylogeny 
	[CAFE assumes that gene birth and death are equally probable, see Hahn et al. (2005)].
	If mu is provided, lambda and mu represent the probability of gene birth and gene death, respectively.
**/
struct square_matrix* compute_birthdeath_rates(double branchlength, double lambda, double mu, int maxFamilysize)
{
	struct square_matrix* matrix = (struct square_matrix*)memory_new(1,sizeof(struct square_matrix)  );
	int sz = maxFamilysize + 1;
	square_matrix_init(matrix, sz);

	square_matrix_set(matrix, 0, 0, 1);		//Once you are zero you are almost surely zero

	double alpha = 0;
	double beta = 0;
	double coeff = 1;

	if (mu < 0 || lambda == mu) {
		alpha = lambda*branchlength/(1+lambda*branchlength);
		beta = alpha;
		coeff = 1 - 2 * alpha;
	}		
	else {
		double e_diff = exp((lambda - mu)*branchlength);
		double numerator = e_diff - 1;
		double denominator = lambda*(e_diff) - mu;
		alpha = (mu*numerator)/denominator;
		beta = (lambda*numerator)/denominator;
		coeff = 1 - alpha - beta;
	}

  init_matrix(matrix, coeff);

	if (coeff > 0 && coeff != 1)
	{
		for (int s = 1 ; s <= maxFamilysize; s++ )
		{ 
			for (int c = 0 ; c <= maxFamilysize ; c++ )
			{
				if (mu < 0)
					square_matrix_set(matrix, s, c, birthdeath_rate_with_log_alpha(s, c, log(alpha), coeff, &cache));
				else
					square_matrix_set(matrix, s, c, birthdeath_rate_with_log_alpha_beta(s, c, log(alpha), log(beta), log(coeff), &cache));
			}
		}
	}
#ifdef VERBOSE
    else
    {
        fprintf(stderr, "WARNING: Zero matrix set for branch %f, lambda %f, mu %f\n", branchlength, lambda, mu);
    }
#endif
	return matrix;
}

void birthdeath_cache_matrix_resize(struct square_matrix* matrix, int remaxFamilysize, double branchlength, double lambda, double mu)
{
	double alpha = lambda*branchlength/(1+lambda*branchlength);
	double coeff = 1 - 2 * alpha;
	int old = matrix->size;
	alpha = log(alpha);
	square_matrix_resize(matrix, remaxFamilysize + 1);

	for (int c = old + 1 ; c <= remaxFamilysize ; c++ )
	{
		square_matrix_set(matrix, 0, c, 0);
	}
	for (int s = 1 ; s <= remaxFamilysize; s++ )
	{
		for (int c = old + 1 ; c <= remaxFamilysize ; c++ )
		{
			square_matrix_set(matrix, s, c, birthdeath_rate_with_log_alpha(s,c,alpha,coeff, &cache));
		}
	}
	for (int s = old  + 1 ; s <= remaxFamilysize; s++ )
	{ 
		for (int c = 0; c <= remaxFamilysize ; c++ )
		{
			square_matrix_set(matrix, s, c, birthdeath_rate_with_log_alpha(s, c, alpha, coeff, &cache));
		}
	}
}

void birthdeath_cache_resize(pBirthDeathCacheArray pbdc_array, int remaxFamilysize)
{
	if (pbdc_array->maxFamilysize >= remaxFamilysize) return;
	chooseln_cache_resize2(&cache, remaxFamilysize);
	void** keys = NULL;
	int num = (int)hash_table_get_keys(pbdc_array->table, &keys);
	for (int i = 0; i<num; i++) {
		struct square_matrix* matrix = hash_table_lookup(pbdc_array->table, keys[i], sizeof(struct BirthDeathCacheKey));
		struct BirthDeathCacheKey* key = (struct BirthDeathCacheKey*)keys[i];
		if (matrix == NULL) continue;
		birthdeath_cache_matrix_resize(matrix, remaxFamilysize, key->branchlength, key->lambda, key->mu);
	}
	pbdc_array->maxFamilysize = remaxFamilysize;
}

pBirthDeathCacheArray birthdeath_cache_init(int size)
{
	pBirthDeathCacheArray pbdc_array = (pBirthDeathCacheArray)memory_new(1, sizeof(BirthDeathCacheArray));
	pbdc_array->table = hash_table_new(MODE_VALUEREF);
	pbdc_array->maxFamilysize = size;

	if (!chooseln_is_init2(&cache))
		chooseln_cache_init2(&cache, pbdc_array->maxFamilysize);
	else if (cache.size < pbdc_array->maxFamilysize)
		chooseln_cache_resize2(&cache, pbdc_array->maxFamilysize);

	return pbdc_array;
}

void birthdeath_cache_array_free(pBirthDeathCacheArray pbdc_array)
{
    void** keys = NULL;
    int num = (int)hash_table_get_keys(pbdc_array->table, &keys);
    for (int i=0; i<num; i++) {
		struct square_matrix* matrix = hash_table_lookup(pbdc_array->table, keys[i], sizeof(struct BirthDeathCacheKey));
		square_matrix_delete(matrix);
    }
    free(keys);
    hash_table_delete(pbdc_array->table);
	memory_free(pbdc_array);
}


/** 
	Returns square matrix of doubles, rows and columns representing the transition probability in birthdeath rate
	in a change of one family size to another, with the given values of branch length, lambda, and mu
**/
struct square_matrix* birthdeath_cache_get_matrix(pBirthDeathCacheArray pbdc_array, double branchlength, double lambda, double mu )
{
	struct BirthDeathCacheKey key;
    memset(&key, 0, sizeof(struct BirthDeathCacheKey));
	key.branchlength = branchlength;
	key.lambda = lambda;
	key.mu = mu;

	struct square_matrix* matrix = hash_table_lookup(pbdc_array->table, &key, sizeof(struct BirthDeathCacheKey));
	if (matrix == NULL)
	{
#ifdef VERBOSE
    if (lambda < 0.000000000001)
      printf("WARNING: building matrix for 0 lambda\n");
#endif
    matrix = compute_birthdeath_rates(key.branchlength, key.lambda, key.mu, pbdc_array->maxFamilysize);
		hash_table_add(pbdc_array->table, &key, sizeof(struct BirthDeathCacheKey), matrix, sizeof(struct square_matrix*));
	}
	return matrix;
}

