#include <assert.h>

#include "birthdeath.h"
#include<utils.h>
#include<stdlib.h>
#include<stdio.h>
#include<pthread.h>
#include<family.h>
#include<float.h>
#include<mathfunc.h>
#include "chooseln_cache.h"

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
	matrix->values = (double**)memory_new_2dim(sz, sz, sizeof(double));
	matrix->size = sz;
}

void square_matrix_set(struct square_matrix* matrix, int x, int y, double val)
{
	assert(x < matrix->size);
	assert(y < matrix->size);
	matrix->values[x][y] = val;
}

void square_matrix_delete(struct square_matrix* matrix)
{
	memory_free_2dim((void**)matrix->values, matrix->size, matrix->size, NULL);
}

void square_matrix_resize(struct square_matrix* matrix, int new_size)
{
	matrix->values = (double**)memory_realloc(matrix->values, new_size, sizeof(double*));

	for (int s = 0; s < matrix->size; s++)
	{
		matrix->values[s] = (double*)memory_realloc(matrix->values[s], new_size, sizeof(double));
	}
	for (int s = matrix->size; s < new_size; s++)
	{
		matrix->values[s] = (double*)memory_new(new_size, sizeof(double));
	}
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
pBirthDeathCache birthdeath_cache_new(double branchlength, double lambda, double mu, int maxFamilysize)
{
	pBirthDeathCache pbdc = (pBirthDeathCache)memory_new(1,sizeof(BirthDeathCache)  );
	int sz = maxFamilysize + 1;
	square_matrix_init(&pbdc->matrix, sz);
	pbdc->branchlength = branchlength;
	pbdc->maxFamilysize = maxFamilysize;
	pbdc->lambda = lambda;
	pbdc->mu = mu;

	square_matrix_set(&pbdc->matrix, 0, 0, 1);		//Once you are zero you are almost surely zero

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

	init_matrix(&pbdc->matrix, coeff);

	if (coeff > 0 && coeff != 1)
	{
		for (int s = 1 ; s <= maxFamilysize; s++ )
		{ 
			for (int c = 0 ; c <= maxFamilysize ; c++ )
			{
				if (mu < 0)
					square_matrix_set(&pbdc->matrix, s, c, birthdeath_rate_with_log_alpha(s, c, log(alpha), coeff, &cache));
				else
					square_matrix_set(&pbdc->matrix, s, c, birthdeath_rate_with_log_alpha_beta(s, c, log(alpha), log(beta), log(coeff), &cache));
			}
		}
	}
	return pbdc;
}


pBirthDeathCache eq_birthdeath_cache_new(double branchlength, double lambda, int maxFamilysize)
{
	return birthdeath_cache_new(branchlength, lambda, -1, maxFamilysize);
}


pBirthDeathCache birthdeath_cach_resize(pBirthDeathCache pbdc, int remaxFamilysize)
{
	int s,c;
	double lambda = pbdc->lambda;
	double branchlength = pbdc->branchlength;

	double alpha = lambda*branchlength/(1+lambda*branchlength);
	double coeff = 1 - 2 * alpha;
	int old = pbdc->maxFamilysize;
	alpha = log(alpha);
	square_matrix_resize(&pbdc->matrix, remaxFamilysize + 1);

	for ( c = old + 1 ; c <= remaxFamilysize ; c++ )
	{
		square_matrix_set(&pbdc->matrix, 0, c, 0);
	}
	for ( s = 1 ; s <= remaxFamilysize; s++ )
	{
		for ( c = old + 1 ; c <= remaxFamilysize ; c++ )
		{
			square_matrix_set(&pbdc->matrix, s, c, birthdeath_rate_with_log_alpha(s,c,alpha,coeff, &cache));
		}
	}
	for ( s = old  + 1 ; s <= remaxFamilysize; s++ )
	{ 
		for ( c = 0; c <= remaxFamilysize ; c++ )
		{
			square_matrix_set(&pbdc->matrix, s, c, birthdeath_rate_with_log_alpha(s, c, alpha, coeff, &cache));
		}
	}
	pbdc->maxFamilysize = remaxFamilysize;
	return pbdc;
}


void birthdeath_cache_free(void* ptr)
{
	pBirthDeathCache pbdc = (pBirthDeathCache)ptr;
	square_matrix_delete(&pbdc->matrix);
    pbdc->matrix.values = NULL;
	memory_free(pbdc);	
	pbdc = NULL;
}

pBirthDeathCache eq_birthdeath_search_list_for_lambda(pArrayList plist, double lambda)
{
	int i;
	for ( i = 0 ; i < plist->size ; i++ )
	{
		pBirthDeathCache pbdc = (pBirthDeathCache)plist->array[i];
		if ( pbdc->lambda == lambda ) return pbdc;
	}
	return NULL;
}

pBirthDeathCache birthdeath_search_list_for_lambda_mu(pArrayList plist, double lambda, double mu)
{
	if (mu < 0) {
		return eq_birthdeath_search_list_for_lambda(plist, lambda);
	}
	int i;
	for ( i = 0 ; i < plist->size ; i++ )
	{
		pBirthDeathCache pbdc = (pBirthDeathCache)plist->array[i];
		if ( pbdc->lambda == lambda && pbdc->mu == mu) {
			return pbdc;
		}
	}
	return NULL;
}


typedef struct
{
	double branchlength;
	double lambda;
	double mu;
	int maxFamilysize;
	pBirthDeathCache pbdc;
}BDCThread;
typedef BDCThread* pBDCThread;

void* __cafe_set_birthdeath_cache_thread_func(void* ptr)
{
	pBDCThread pbdt = (pBDCThread)ptr;
	pbdt->pbdc = birthdeath_cache_new(pbdt->branchlength, pbdt->lambda, pbdt->mu, pbdt->maxFamilysize);
	return (NULL);
}

void cafe_set_birthdeath_cache_thread(pCafeTree tree, int k_value, int* family_sizes, int* rootfamily_sizes)
{
	if ( tree->pbdc_array )
	{
		birthdeath_cache_array_free( tree->pbdc_array );
	}
	pBirthDeathCacheArray pbdc_array = (pBirthDeathCacheArray)memory_new(1,sizeof(BirthDeathCacheArray));
	pbdc_array->maxFamilysize = MAX(family_sizes[1], rootfamily_sizes[1]);
	
	if ( !chooseln_is_init() ) 
	{
		chooseln_cache_init2(&cache, pbdc_array->maxFamilysize );
	}
	else if ( cache.size < pbdc_array->maxFamilysize ) 
	{
		chooseln_cache_resize2(&cache, pbdc_array->maxFamilysize );
	}
	
	int i,j,k,l = 0;
	pArrayList nlist = ((pTree)tree)->nlist;
	pArrayList thread_param = arraylist_new( nlist->size+1 );	// thread_param is an array of void* size of nodes(+1). 
	
	for( i = 0, j = 0 ; j < nlist->size; j++ )
	{
		pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[j];
		pCafeNode pcnode = (pCafeNode)pnode;
		if ( pnode->branchlength <= 0 ) continue;
		//param->branchlengths_sorted[i++] = pnode->branchlength;
		
		// look for param_lambdas and param_mus 
		if ( pcnode->param_lambdas ) {
			if (pcnode->param_mus) {
				for ( k=0; k < k_value; k++) {
					for ( l = 0 ; l < thread_param->size ; l++ )				// at first thread_param->size is zero, doesn't go into this loop
					{
						pBDCThread pbdt = (pBDCThread)thread_param->array[l];	// previously added thread
						if ( pbdt->branchlength == pnode->branchlength 
							 && pbdt->lambda == pcnode->param_lambdas[k]
							 && pbdt->mu == pcnode->param_mus[k] )	
						{
							break;		// if there exist a thread with branchlength and lambda of pnode already
						}
					}
					if ( l == thread_param->size )	// if thread of pnode is not found make thread of pnode
					{
						pBDCThread pbdt = (pBDCThread)memory_new(1,sizeof(BDCThread));
						pbdt->branchlength = pnode->branchlength;
						pbdt->lambda = pcnode->param_lambdas[k];
						pbdt->mu = pcnode->param_mus[k];
						pbdt->maxFamilysize = pbdc_array->maxFamilysize;
						pbdt->pbdc = NULL;
						arraylist_add( thread_param, pbdt ); // add new thread to thread_param
					}
				}
			}
			else {
				for ( k=0; k < k_value; k++) {
					for ( l = 0 ; l < thread_param->size ; l++ )				// at first thread_param->size is zero, doesn't go into this loop
					{
						pBDCThread pbdt = (pBDCThread)thread_param->array[l];	// previously added thread
						if ( pbdt->branchlength == pnode->branchlength 
							 && pbdt->lambda == pcnode->param_lambdas[k]
							 && pbdt->mu == pcnode->mu)
						{
							break;		// if there exist a thread with branchlength and lambda of pnode already
						}
					}
					if ( l == thread_param->size )	// if thread of pnode is not found make thread of pnode
					{
						pBDCThread pbdt = (pBDCThread)memory_new(1,sizeof(BDCThread));
						pbdt->branchlength = pnode->branchlength;
						pbdt->lambda = pcnode->param_lambdas[k];
						pbdt->mu = pcnode->mu;
						pbdt->maxFamilysize = pbdc_array->maxFamilysize;
						pbdt->pbdc = NULL;
						arraylist_add( thread_param, pbdt ); // add new thread to thread_param
					}
				}
			}
		}
    else 
    {
      // look for thread with branchlength and lambda and mu of pnode in the array of thread_param
      for ( k = 0 ; k < thread_param->size ; k++ )				// at first thread_param->size is zero, doesn't go into this loop
      {
        pBDCThread pbdt = (pBDCThread)thread_param->array[k];	// previously added thread
        if ( pbdt->branchlength == pnode->branchlength 
             && pbdt->lambda == pcnode->lambda 
             && pbdt->mu == pcnode->mu )	
        {
          break;		// if there exist a thread with branchlength and lambda of pnode already
        }
      }
      if ( k == thread_param->size )	// if thread of pnode is not found make thread of pnode
      {
        pBDCThread pbdt = (pBDCThread)memory_new(1,sizeof(BDCThread));
        pbdt->branchlength = pnode->branchlength;
        pbdt->lambda = ((pCafeNode)pnode)->lambda;
        pbdt->mu = ((pCafeNode)pnode)->mu;
        pbdt->maxFamilysize = pbdc_array->maxFamilysize;
        pbdt->pbdc = NULL;
        arraylist_add( thread_param, pbdt ); // add new thread to thread_param
      }
    }
		
	}
	// now there is thread for every node
	
	//qsort(param->branchlengths_sorted, param->num_branches, sizeof(int), __cmp_int );
	int numthreads = thread_param->size;
	thread_run_with_arraylist(numthreads, __cafe_set_birthdeath_cache_thread_func, thread_param ); // run function on array of threads
    pbdc_array->table = hash_table_new(MODE_VALUEREF);
    for ( i = 0 ; i < thread_param->size ; i++ )                // for each thread existing for all (branch lengths, parameters) combinations
	{
		pBDCThread pbdt = (pBDCThread)thread_param->array[i];
        double* key = &pbdt->branchlength;
		pArrayList plist = (pArrayList) hash_table_lookup(pbdc_array->table, key, sizeof(double));
        if ( plist == NULL )
		{
            plist = arraylist_new(10);
			arraylist_add( plist, pbdt->pbdc );
			hash_table_add(pbdc_array->table, (void*)key, sizeof(double), (void *)plist, sizeof(ArrayList));
			
		}
		else
		{
			if ( birthdeath_search_list_for_lambda_mu(plist,pbdt->lambda, pbdt->mu) == NULL )
			{
				arraylist_add( plist, pbdt->pbdc );
			}
			
		}
	}
	
	arraylist_free(thread_param,free);
	tree->pbdc_array = pbdc_array;
	cafe_tree_set_birthdeath(tree);
}

void cafe_set_birthdeath_cache(pCafeParam param)
/* 
 * Before: param->param_set_func(param, param->lambda);
 */
{
	if ( param->pcafe->pbdc_array )
	{
		birthdeath_cache_array_free( param->pcafe->pbdc_array );
	}

	pBirthDeathCacheArray pbdc_array = (pBirthDeathCacheArray)memory_new(1,sizeof(BirthDeathCacheArray));

	int i;
	pArrayList nlist = ((pTree)param->pcafe)->nlist;

    pbdc_array->table = hash_table_new(MODE_VALUEREF);
	pbdc_array->maxFamilysize = MAX(param->family_sizes[1], param->rootfamily_sizes[1]);

	if ( !chooseln_is_init2(&cache) ) 
		chooseln_cache_init2(&cache, pbdc_array->maxFamilysize );
	else if ( cache.size < pbdc_array->maxFamilysize ) 
		chooseln_cache_resize2(&cache, pbdc_array->maxFamilysize );

	for ( i = 0 ; i < nlist->size ; i++ )
	{
		pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
		if ( pnode->branchlength > 0 ) 
		{	
            double* key = &pnode->branchlength;
			pArrayList plist = hash_table_lookup(pbdc_array->table, key, sizeof(double));
			double lambda = ((pCafeNode)pnode)->lambda;
			double mu = ((pCafeNode)pnode)->mu;
			if (plist == NULL)
            {
				plist = arraylist_new(10);
				pBirthDeathCache cache = birthdeath_cache_new(pnode->branchlength, lambda, mu, pbdc_array->maxFamilysize);
				arraylist_add( plist, cache);
                hash_table_add(pbdc_array->table, key, sizeof(double), plist, sizeof(ArrayList));
			}
			else
			{
				if ( birthdeath_search_list_for_lambda_mu(plist, lambda, mu) == NULL )
				{
					pBirthDeathCache cache = birthdeath_cache_new(pnode->branchlength, lambda, mu, pbdc_array->maxFamilysize);
					arraylist_add( plist, cache);
				}
			}
		}
	}
	param->pcafe->pbdc_array = pbdc_array;
	cafe_tree_set_birthdeath(param->pcafe);
}

void cafe_resize_birthdeath_cache(pCafeParam param)
{
	pBirthDeathCacheArray pbdc_array = param->pcafe->pbdc_array;
	int remaxFamilysize = MAX(param->family_sizes[1], param->rootfamily_sizes[1]);
	if ( pbdc_array->maxFamilysize >= remaxFamilysize ) return;
	chooseln_cache_resize2(&cache, remaxFamilysize);
	int i,j;
    void** keys = NULL;
    int num = (int)hash_table_get_keys(pbdc_array->table, &keys);
    for (i=0; i<num; i++) {
        pArrayList plist = hash_table_lookup(pbdc_array->table, keys[i], sizeof(double));
		if ( plist == NULL ) continue;
		for ( j = 0 ; j < plist->size ; j++ )
		{
			pBirthDeathCache pbdc = (pBirthDeathCache)plist->array[j];
			birthdeath_cach_resize(pbdc, remaxFamilysize);
		}
    }
	pbdc_array->maxFamilysize = remaxFamilysize;
	cafe_tree_set_birthdeath(param->pcafe);
}

void birthdeath_cache_array_free(pBirthDeathCacheArray pbdc_array)
{
	int i ;
    void** keys = NULL;
    int num = (int)hash_table_get_keys(pbdc_array->table, &keys);
    for (i=0; i<num; i++) {
        pArrayList plist = hash_table_lookup(pbdc_array->table, keys[i], sizeof(double));
        arraylist_free( plist, birthdeath_cache_free );
    }
    hash_table_delete(pbdc_array->table);
	memory_free(pbdc_array);
	pbdc_array = NULL;
}


struct square_matrix* eq_birthdeath_cache_get_matrix(pBirthDeathCacheArray pbdc_array, double branchlength, double lambda )
{
	pArrayList plist;
	pBirthDeathCache pbdc = NULL;
    double* key = &branchlength;
    plist = (pArrayList)hash_table_lookup(pbdc_array->table, key, sizeof(double));
	if ( plist == NULL )
	{
		plist = (pArrayList)arraylist_new(10);
		pbdc = eq_birthdeath_cache_new( branchlength , lambda, pbdc_array->maxFamilysize );
		arraylist_add(plist, pbdc);
        hash_table_add(pbdc_array->table, key, sizeof(double), plist, sizeof(ArrayList));
	}
	else if ( (pbdc = eq_birthdeath_search_list_for_lambda(plist,lambda)) == NULL )
	{
		pbdc = eq_birthdeath_cache_new(branchlength, lambda, pbdc_array->maxFamilysize );
		arraylist_add( plist, pbdc );
	}
	if ( pbdc == NULL )
	{
		plist = (pArrayList)arraylist_new(10);
		pbdc = eq_birthdeath_cache_new( branchlength , lambda, pbdc_array->maxFamilysize );
		arraylist_add(plist, pbdc);
        hash_table_add(pbdc_array->table, key, sizeof(double), plist, sizeof(ArrayList));
	}
	return &pbdc->matrix;
}

/** 
	Returns square matrix of doubles, rows and columns representing the transition probability in birthdeath rate
	in a change of one family size to another, with the given values of branch length, lambda, and mu
**/
struct square_matrix* birthdeath_cache_get_matrix(pBirthDeathCacheArray pbdc_array, double branchlength, double lambda, double mu )
{
	if (mu < 0) {
		return eq_birthdeath_cache_get_matrix(pbdc_array, branchlength, lambda);
	}
	pArrayList plist;
	pBirthDeathCache pbdc = NULL;
    double* key = &branchlength;
    plist = (pArrayList) hash_table_lookup(pbdc_array->table, key, sizeof(double));
	if ( plist == NULL )
	{
		plist = (pArrayList)arraylist_new(10);
		pbdc = birthdeath_cache_new( branchlength, lambda, mu, pbdc_array->maxFamilysize );
		arraylist_add(plist, pbdc);
        hash_table_add(pbdc_array->table, key, sizeof(double), plist, sizeof(ArrayList));
	}
	else if ( (pbdc = birthdeath_search_list_for_lambda_mu(plist,lambda, mu)) == NULL )
	{
		pbdc = birthdeath_cache_new(branchlength, lambda, mu, pbdc_array->maxFamilysize );
		arraylist_add( plist, pbdc );
	}
	if ( pbdc == NULL )
	{
		plist = (pArrayList)arraylist_new(10);
		pbdc = birthdeath_cache_new( branchlength, lambda, mu, pbdc_array->maxFamilysize );
		arraylist_add(plist, pbdc);
        hash_table_add(pbdc_array->table, key, sizeof(double), plist, sizeof(ArrayList));
	}
	return &pbdc->matrix;
}

