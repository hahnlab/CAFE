#include "birthdeath.h"
#include<utils.h>
#include<stdlib.h>
#include<stdio.h>
#include<pthread.h>
#include<family.h>
#include<float.h>
#include<mathfunc.h>

/*
P(X(t) = c | X(0) = s)  = \sum_{j=0}^{\min(s,c)} \binom{s}{j}\binom{s+c-j-1}{s-1}
                                                 \alpha^{s+c-2j}(1-2\alpha)^j
 */

int chooseln_cache_size;
static double** chooseln_cache;

int chooseln_is_init()
{
	return chooseln_cache ? 1 : 0;
}

int get_chooseln_cache_size() 
{ 
	return chooseln_cache_size; 
}

double chooseln_get(int n, int x)
{
	if ( chooseln_cache[n] && chooseln_cache[n][x] >= 0 ) return chooseln_cache[n][x];
	if ( chooseln_cache[n] == NULL )
	{
		int i;
		chooseln_cache[n] = (double*)memory_new(chooseln_cache_size+1, sizeof(double) );
		for( i = 0 ; i <= chooseln_cache_size ; i++ ) chooseln_cache[n][i] = -1.0;	
	}
	chooseln_cache[n][x] = chooseln(n,x);	
	return chooseln_cache[n][x];
}

void chooseln_cache_preset(int maxFamilysize, int sfrom )
{
	int s, c, j;
	for ( s = sfrom; s <= maxFamilysize; s++ )
	{ 
		for ( c = 0; c <= maxFamilysize ; c++ )
		{
			int m = MIN(s,c);
			double s_add_c = s + c;
			double s_add_c_sub_1 = s_add_c - 1;
			double s_sub_1 = s - 1;
			for ( j = 0 ; j <= m ; j++ )
			{
				chooseln_get(s,j); 
				chooseln_get(s_add_c_sub_1-j,s_sub_1);
			}
		}
	}
}

void chooseln_cache_resize(int resize)
{
	if ( chooseln_cache_size >= resize ) return;
	int i;
	if ( chooseln_cache )
	{
		chooseln_cache = (double**)memory_realloc(chooseln_cache,resize*2,sizeof(double*));
	}
	else
	{
		chooseln_cache = (double**)memory_new(resize*2,sizeof(double*));
	}
	for ( i = chooseln_cache_size*2 ; i < resize*2 ; i++ )
	{
		chooseln_cache[i] = NULL;
	}
	int oldsize = chooseln_cache_size;
	chooseln_cache_size = resize;
	chooseln_cache_preset( resize, oldsize+1);
	fprintf(stderr, "** Cache resize: %d ==> %d\n", oldsize, resize );
}

void chooseln_cache_init(int size)
{
	chooseln_cache_size = size;
	chooseln_cache = (double**)memory_new(size*2,sizeof(double*));
	int i;
	for( i = 0; i < size*2 ; i++ ) 
	{
		chooseln_cache[i] = NULL;
	}
	chooseln_cache_preset(size, 1);
}

void chooseln_cache_free()
{
	int i;
	for( i = 0; i < chooseln_cache_size*2 ; i++ )
	{
		if ( chooseln_cache[i] ) memory_free(chooseln_cache[i]);
		chooseln_cache[i] = NULL;
	}
	memory_free(chooseln_cache);
	chooseln_cache = NULL;
}




double birthdeath_rate_with_log_alpha_beta(int s, int c, double log_alpha, double log_beta, double log_coeff )
{
	int j;
	int m = MIN(c,s);
	double t, p = 0;
	int s_add_c = s + c;
	int s_add_c_sub_1 = s_add_c - 1;
	int s_sub_1 = s - 1;
	for ( j = 0, p = 0 ; j <= m ; j++ )
	{
		t = chooseln_cache[s][j] + chooseln_cache[s_add_c_sub_1-j][s_sub_1] + (s-j)*log_alpha + (c-j)*log_beta + j*log_coeff;
		//t = chooseln_get(s, j) + chooseln_get(s_add_c_sub_1-j,s_sub_1) + (s-j)*log_alpha + (c-j)*log_beta + j*log_coeff;
		p += exp(t);
	}
	return MAX(MIN(p,1),0);
}

double birthdeath_rate_with_log_alpha(int s, int c, double log_alpha, double coeff )
{
	int j;
	int m = MIN(c,s);
	double t, p, lastterm = 1;
	int s_add_c = s + c;
	int s_add_c_sub_1 = s_add_c - 1;
	int s_sub_1 = s - 1;
	for ( j = 0, p = 0 ; j <= m ; j++ )
	{
		t = chooseln_cache[s][j] + chooseln_cache[s_add_c_sub_1-j][s_sub_1] + (s_add_c-2*j)*log_alpha;
		//t = chooseln_get(s, j) + chooseln_get(s_add_c_sub_1-j,s_sub_1) + (s_add_c-2*j)*alpha;
		p += (exp(t) * lastterm);
		lastterm *= coeff;
	}
	return MAX(MIN(p,1),0);
}

double birthdeath_likelihood_with_s_c(int s, int c, double branchlength, double lambda, double mu) 
{	
	double alpha, coeff, beta=0;
	double denominator, numerator = 0;	
	
	if (s == 0) {
		if (c == 0) {
			return 1;
		}
		else {
			return 0;
		}
	}
	else {
		if ((mu < 0) || (lambda == mu)) {
			alpha = lambda*branchlength/(1+lambda*branchlength);
			coeff = 1 - 2 * alpha;
			if (coeff <= 0) {
				return 0;
			}
			else {
				return birthdeath_rate_with_log_alpha(s,c,log(alpha),coeff);
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
				return birthdeath_rate_with_log_alpha_beta(s,c,log(alpha),log(beta),log(coeff));
			}
		}
	}
}

// THE FUNCTION!!!!!!!!!!
// must add mu to calculate the transition probability 
pBirthDeathCache birthdeath_cache_new(double branchlength, double lambda, double mu, int maxFamilysize )
{
	int s, c;
	double alpha = 0;
	double beta = 0;
	double coeff = 1;
	double denominator = 0;
	double numerator = 0;

	if (mu < 0) {
		return eq_birthdeath_cache_new(branchlength, lambda, maxFamilysize);
	}
	pBirthDeathCache pbdc = (pBirthDeathCache)memory_new(1,sizeof(BirthDeathCache)  );
	pbdc->matrix = (double**)memory_new_2dim( (maxFamilysize+1) , (maxFamilysize+1), sizeof(double) );
	pbdc->attrib.branchlength = branchlength;
	pbdc->attrib.maxFamilysize = maxFamilysize;
	pbdc->attrib.lambda = lambda;
	pbdc->attrib.mu = mu;

#ifdef __DEBUG__
//	fprintf(stderr, "ADD Branch in cache(%d): %d, %f\n", maxFamilysize, branchlength, lambda);
#endif

	pbdc->matrix[0][0] = 1;		//Once you are zero you are almost surely zero

	if (lambda == mu) {
		alpha = lambda*branchlength/(1+lambda*branchlength);
		beta = alpha;
		coeff = 1 - 2 * alpha;
	}		
	else {
		denominator = lambda*(exp((lambda-mu)*branchlength))-mu;
		numerator = exp((lambda-mu)*branchlength)-1;
		alpha = (mu*numerator)/denominator;
		beta = (lambda*numerator)/denominator;
		coeff = 1 - alpha - beta;
	}

	for ( c = 1 ; c <= maxFamilysize; c++ )
	{
		pbdc->matrix[0][c] = 0;
	}
  if ( coeff <= 0 ) 
  {
    for ( s = 1 ; s <= maxFamilysize; s++ )
    { 
      for ( c = 0 ; c <= maxFamilysize ; c++ )
      {
        pbdc->matrix[s][c] = 0;
      }
    }
  }
  else if ( coeff == 1) 
  {   
    for ( s = 1 ; s <= maxFamilysize; s++ )
    { 
      for ( c = 0 ; c <= maxFamilysize ; c++ )
      {
        if (s == c) {
          pbdc->matrix[s][c] = 1;
        }
        else {
          pbdc->matrix[s][c] = 0;
        }
      }
    }
  }
	else
	{
		for ( s = 1 ; s <= maxFamilysize; s++ )
		{ 
			for ( c = 0 ; c <= maxFamilysize ; c++ )
			{
				pbdc->matrix[s][c] = birthdeath_rate_with_log_alpha_beta(s,c,log(alpha), log(beta), log(coeff));
			}
		}
	}
	return pbdc;
}

pBirthDeathCache eq_birthdeath_cache_new(double branchlength, double lambda, int maxFamilysize )
{
	pBirthDeathCache pbdc = (pBirthDeathCache)memory_new(1,sizeof(BirthDeathCache)  );
	pbdc->matrix = (double**)memory_new_2dim( (maxFamilysize+1) , (maxFamilysize+1), sizeof(double) );
	pbdc->attrib.branchlength = branchlength;
	pbdc->attrib.maxFamilysize = maxFamilysize;
	pbdc->attrib.lambda = lambda;
	pbdc->attrib.mu = -1;

#ifdef __DEBUG__
//	fprintf(stderr, "ADD Branch in cache(%d): %d, %f\n", maxFamilysize, branchlength, lambda);
#endif

	int s, c;
	pbdc->matrix[0][0] = 1;

	double alpha = lambda*branchlength/(1+lambda*branchlength);
	double coeff = 1 - 2 * alpha;

	for ( c = 1 ; c <= maxFamilysize; c++ )
	{
		pbdc->matrix[0][c] = 0;
	}
	if ( coeff <= 0 ) 
	{
		for ( s = 1 ; s <= maxFamilysize; s++ )
		{ 
			for ( c = 0 ; c <= maxFamilysize ; c++ )
			{
				pbdc->matrix[s][c] = 0;
			}
		}
	}
  else if ( coeff == 1) 
  {   
    for ( s = 1 ; s <= maxFamilysize; s++ )
    { 
      for ( c = 0 ; c <= maxFamilysize ; c++ )
      {
        if (s == c) {
          pbdc->matrix[s][c] = 1;
        }
        else {
          pbdc->matrix[s][c] = 0;
        }
      }
    }
  }
	else
	{
		alpha = log(alpha);
		for ( s = 1 ; s <= maxFamilysize; s++ )
		{ 
			for ( c = 0 ; c <= maxFamilysize ; c++ )
			{
				pbdc->matrix[s][c] = birthdeath_rate_with_log_alpha(s,c,alpha,coeff);
			}
		}
	}
	return pbdc;
}

pBirthDeathCache birthdeath_cach_resize(pBirthDeathCache pbdc, int remaxFamilysize)
{
	int s,c;
	double lambda = pbdc->attrib.lambda;
	double branchlength = pbdc->attrib.branchlength;
#ifdef __DEBUG__
	fprintf(stderr, "Increae Branch family size from %d -> %d in cache %d, %f\n", pbdc->attrib.maxFamilysize, remaxFamilysize, branchlength, lambda);
#endif
	double alpha = lambda*branchlength/(1+lambda*branchlength);
	double coeff = 1 - 2 * alpha;
	int old = pbdc->attrib.maxFamilysize;
	alpha = log(alpha);
	pbdc->matrix = (double**)memory_realloc( pbdc->matrix, remaxFamilysize + 1, sizeof(double*) );

	for ( s = 0 ; s <= old ; s++ )
	{
		pbdc->matrix[s] = (double*)memory_realloc( pbdc->matrix[s], remaxFamilysize + 1, sizeof(double));
	}
	for ( s = old + 1 ; s <= remaxFamilysize; s++ )
	{
		pbdc->matrix[s] = (double*)memory_new( remaxFamilysize + 1, sizeof(double));
	}

	for ( c = old + 1 ; c <= remaxFamilysize ; c++ )
	{
		pbdc->matrix[0][c] = 0;
	}
	for ( s = 1 ; s <= remaxFamilysize; s++ )
	{
		for ( c = old + 1 ; c <= remaxFamilysize ; c++ )
		{
			pbdc->matrix[s][c] = birthdeath_rate_with_log_alpha(s,c,alpha,coeff);
		}
	}
	for ( s = old  + 1 ; s <= remaxFamilysize; s++ )
	{ 
		for ( c = 0; c <= remaxFamilysize ; c++ )
		{
			pbdc->matrix[s][c] = birthdeath_rate_with_log_alpha(s,c,alpha,coeff);
		}
	}
	pbdc->attrib.maxFamilysize = remaxFamilysize;
	return pbdc;
}


void birthdeath_cache_free(void* ptr)
{
	pBirthDeathCache pbdc = (pBirthDeathCache)ptr;
	memory_free_2dim( (void**)pbdc->matrix, 
			          pbdc->attrib.maxFamilysize + 1, 
					  pbdc->attrib.maxFamilysize + 1, NULL );
    pbdc->matrix = NULL;
	memory_free(pbdc);	
	pbdc = NULL;
}

pBirthDeathCache eq_birthdeath_search_list_for_lambda(pArrayList plist, double lambda)
{
	int i;
	for ( i = 0 ; i < plist->size ; i++ )
	{
		pBirthDeathCache pbdc = (pBirthDeathCache)plist->array[i];
		if ( pbdc->attrib.lambda == lambda ) return pbdc;
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
		if ( pbdc->attrib.lambda == lambda && pbdc->attrib.mu == mu) {
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

/* not needed 
void* __cafe_set_eq_birthdeath_cache_thread_func(void* ptr)
{
	pBDCThread pbdt = (pBDCThread)ptr;
	pbdt->pbdc = eq_birthdeath_cache_new(pbdt->branchlength, pbdt->lambda, pbdt->maxFamilysize);
	return (NULL);
}
*/

void cafe_set_birthdeath_cache_thread(pCafeParam param)
{
	if ( param->pcafe->pbdc_array )
	{
		birthdeath_cache_array_free( param->pcafe->pbdc_array );
	}
	pBirthDeathCacheArray pbdc_array = (pBirthDeathCacheArray)memory_new(1,sizeof(BirthDeathCacheArray));
	pbdc_array->maxFamilysize = MAX(param->family_sizes[1], param->rootfamily_sizes[1]);
	
	if ( !chooseln_is_init() ) 
	{
		chooseln_cache_init( pbdc_array->maxFamilysize );
	}
	else if ( chooseln_cache_size < pbdc_array->maxFamilysize ) 
	{
		chooseln_cache_resize( pbdc_array->maxFamilysize );
	}
	
	int i,j,k,l = 0;
	pArrayList nlist = ((pTree)param->pcafe)->nlist;
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
				for ( k=0; k < param->k; k++) {
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
				for ( k=0; k < param->k; k++) {
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
	//int* bl = param->branchlengths_sorted;
	//int size = param->num_branches;
	//pbdc_array->base_bl = bl[0];
	// here is where to fix integer branch length.
	/*pbdc_array->list = arraylist_new(bl[size-1] - bl[0] + 1 + 100);
    for ( i = bl[0]; i <= bl[size-1]; i++ )
	{
		arraylist_add(pbdc_array->list, NULL );
	}
    for ( i = 0 ; i < thread_param->size ; i++ )
    {
		pBDCThread pbdt = (pBDCThread)thread_param->array[i];
		pArrayList plist;
		int idx = pbdt->branchlength - bl[0];
		if ( (plist=pbdc_array->list->array[idx]) == NULL )
		{
			plist = arraylist_new(10);
			arraylist_add( plist, pbdt->pbdc );
			pbdc_array->list->array[idx] = plist;
		}
		else
		{
			if ( birthdeath_search_list_for_lambda_mu(plist,pbdt->lambda, pbdt->mu) == NULL )
			{
				arraylist_add( plist, pbdt->pbdc );
			}
			
		}
	}*/
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
	param->pcafe->pbdc_array = pbdc_array;
	cafe_tree_set_birthdeath(param->pcafe);
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

	/*for( i = 0, j = 0 ; j < nlist->size; j++ )
	{
		pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[j];
		if ( pnode->branchlength > 0 )
		{
			param->branchlengths_sorted[i++] = pnode->branchlength;
		}
	}
	qsort(param->branchlengths_sorted, param->num_branches, sizeof(int), __cmp_int );
     */
	//int* bl = param->branchlengths_sorted;
	//int size = param->num_branches;
	//pbdc_array->base_bl = bl[0];
	//pbdc_array->list = arraylist_new(bl[size-1] - bl[0] + 1 + 100);
    pbdc_array->table = hash_table_new(MODE_VALUEREF);
	pbdc_array->maxFamilysize = MAX(param->family_sizes[1], param->rootfamily_sizes[1]);

	if ( !chooseln_is_init() ) chooseln_cache_init( pbdc_array->maxFamilysize );
	else if ( chooseln_cache_size < pbdc_array->maxFamilysize ) chooseln_cache_resize( pbdc_array->maxFamilysize );

	/*for ( i = bl[0]; i <= bl[size-1]; i++ )
	{
		arraylist_add(pbdc_array->list, NULL );
	}*/

	for ( i = 0 ; i < nlist->size ; i++ )
	{
		pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
		if ( pnode->branchlength > 0 ) 
		{	
            double* key = &pnode->branchlength;
			pArrayList plist = hash_table_lookup(pbdc_array->table, key, sizeof(double));
			double lambda = ((pCafeNode)pnode)->lambda;
			double mu = ((pCafeNode)pnode)->mu;
			//if ( (plist=pbdc_array->list->array[(int)pnode->branchlength-bl[0]]) == NULL )
			if (plist == NULL)
            {
				plist = arraylist_new(10);
				arraylist_add( plist, birthdeath_cache_new(pnode->branchlength, lambda, mu, pbdc_array->maxFamilysize));
				//pbdc_array->list->array[(int)pnode->branchlength-bl[0]] = plist;
                hash_table_add(pbdc_array->table, key, sizeof(double), plist, sizeof(ArrayList));
			}
			else
			{
				if ( birthdeath_search_list_for_lambda_mu(plist, lambda, mu) == NULL )
				{
					arraylist_add( plist, birthdeath_cache_new(pnode->branchlength, lambda, mu, pbdc_array->maxFamilysize));
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
	chooseln_cache_resize(remaxFamilysize);
	int i,j;
	/*for( i = 0 ; i < pbdc_array->list->size ; i++ )
	{
		pArrayList plist = (pArrayList)pbdc_array->list->array[i];
		if ( plist == NULL ) continue;
		for ( j = 0 ; j < plist->size ; j++ )
		{
			pBirthDeathCache pbdc = (pBirthDeathCache)plist->array[j];
			birthdeath_cach_resize(pbdc, remaxFamilysize);
		}
	}*/
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
	/*for ( i = 0 ; i <  pbdc_array->list->size ; i++ )
	{
		if ( pbdc_array->list->array[i] )
		{
			arraylist_free( (pArrayList)pbdc_array->list->array[i], birthdeath_cache_free );
            pbdc_array->list->array[i] = NULL;
		}
	}
	arraylist_free( pbdc_array->list, NULL );
    pbdc_array->list = NULL;*/
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


double** eq_birthdeath_cache_get_matrix(pBirthDeathCacheArray pbdc_array, double branchlength, double lambda )
{
	//int i = branchlength - pbdc_array->base_bl;
	pArrayList plist;
	pBirthDeathCache pbdc = NULL;
	//if ( i < pbdc_array->list->size && i >= 0 )
	//{
		//plist = (pArrayList)pbdc_array->list->array[i];
    double* key = &branchlength;
    plist = (pArrayList)hash_table_lookup(pbdc_array->table, key, sizeof(double));
		if ( plist == NULL )
		{
			plist = (pArrayList)arraylist_new(10);
			pbdc = eq_birthdeath_cache_new( branchlength , lambda, pbdc_array->maxFamilysize );
			arraylist_add(plist, pbdc);
			//pbdc_array->list->array[i] = plist;
            hash_table_add(pbdc_array->table, key, sizeof(double), plist, sizeof(ArrayList));
		}
		else if ( (pbdc = eq_birthdeath_search_list_for_lambda(plist,lambda)) == NULL )
		{
			pbdc = eq_birthdeath_cache_new(branchlength, lambda, pbdc_array->maxFamilysize );
			arraylist_add( plist, pbdc );
		}
	//}
	if ( pbdc == NULL )
	{
		/*for ( i = pbdc_array->base_bl + pbdc_array->list->size ;  i <= branchlength ; i++ )
		{
			arraylist_add(pbdc_array->list,NULL);
		}*/
		plist = (pArrayList)arraylist_new(10);
		pbdc = eq_birthdeath_cache_new( branchlength , lambda, pbdc_array->maxFamilysize );
		arraylist_add(plist, pbdc);
		//pbdc_array->list->array[i - pbdc_array->base_bl - 1] = plist;
        hash_table_add(pbdc_array->table, key, sizeof(double), plist, sizeof(ArrayList));
	}
	return pbdc->matrix;
}

double** birthdeath_cache_get_matrix(pBirthDeathCacheArray pbdc_array, double branchlength, double lambda, double mu )
{
	if (mu < 0) {
		return eq_birthdeath_cache_get_matrix(pbdc_array, branchlength, lambda);
	}
	//int i = branchlength - pbdc_array->base_bl;
	pArrayList plist;
	pBirthDeathCache pbdc = NULL;
	//if ( i < pbdc_array->list->size && i >= 0 )
	//{
    double* key = &branchlength;
    plist = (pArrayList) hash_table_lookup(pbdc_array->table, key, sizeof(double));
		//plist = (pArrayList)pbdc_array->list->array[i];
		if ( plist == NULL )
		{
			plist = (pArrayList)arraylist_new(10);
			pbdc = birthdeath_cache_new( branchlength, lambda, mu, pbdc_array->maxFamilysize );
			arraylist_add(plist, pbdc);
			//pbdc_array->list->array[i] = plist;
            hash_table_add(pbdc_array->table, key, sizeof(double), plist, sizeof(ArrayList));
		}
		else if ( (pbdc = birthdeath_search_list_for_lambda_mu(plist,lambda, mu)) == NULL )
		{
			pbdc = birthdeath_cache_new(branchlength, lambda, mu, pbdc_array->maxFamilysize );
			arraylist_add( plist, pbdc );
		}
	//}
	if ( pbdc == NULL )
	{
		/*for ( i = pbdc_array->base_bl + pbdc_array->list->size ;  i <= branchlength ; i++ )
		{
			arraylist_add(pbdc_array->list,NULL);
		}*/
		plist = (pArrayList)arraylist_new(10);
		pbdc = birthdeath_cache_new( branchlength, lambda, mu, pbdc_array->maxFamilysize );
		arraylist_add(plist, pbdc);
		//pbdc_array->list->array[i - pbdc_array->base_bl - 1] = plist;
        hash_table_add(pbdc_array->table, key, sizeof(double), plist, sizeof(ArrayList));
	}
	return pbdc->matrix;
}

double birthdeath_cache_get(pBirthDeathCacheArray pbdc_array, 
		                    int s, int c, int branchlength, double lambda, double mu )
{
	return birthdeath_cache_get_matrix(pbdc_array,branchlength,lambda,mu)[s][c];
}
