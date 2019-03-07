/*! 
* The variable \ref cafe_param is a global singleton that holds general program state.
* \ref cafe_param holds pcafe, a \ref CafeTree, and pfamily, a \ref CafeFamily . These are set
* by the user via the commands "tree" and "load" respectively. When the user calls
* the command "lambda" calculations are performed on the pcafe and pfamily variables.
*/
#include "cafe.h"
#include<stdlib.h>
#include<math.h>
#include<stdio.h>
#include<pthread.h>
#include<mathfunc.h>
#include<memalloc.h>
#include<utils.h>

pBirthDeathCacheArray probability_cache = NULL;
extern struct chooseln_cache cache;

/**
\brief Logs the message and parameters in a standard way
*
* If the user has assigned a file for logging, message
* is written to the file and also to stdout. The stdout
* write will be suppressed if the param "quiet" flag is set.
*
*/void cafe_log(pCafeParam param, const char* msg, ... )
{
  va_list ap;
  va_start(ap, msg);
  if ( param->flog != stderr && param->flog != stdout )
  {
    va_list ap2; 
    va_copy(ap2, ap); 
    vfprintf(param->flog, msg, ap2); 
    va_end(ap2);
  }
  if (!param->quiet)
  {
	  vfprintf(stdout, msg, ap);
	  fflush(stdout);
  }
  fflush(param->flog);
  va_end(ap);
}

void cafe_free_birthdeath_cache(pCafeTree pcafe)
{
	birthdeath_cache_array_free(probability_cache);
	probability_cache = NULL;
}

void copy_range_to_tree(pCafeTree tree, family_size_range* range)
{
	tree->range.root_min = range->root_min;
	tree->range.root_max = range->root_max;
	tree->range.min = range->min;
	tree->range.max = range->max;
	tree->rfsize = range->root_max - range->root_min + 1;

}

void thread_run(int numthreads, void* (*run)(void*), void* param, int size )
{
	int i;
	pthread_t* pthreads = (pthread_t*)memory_new( numthreads, sizeof(pthread_t));
	for ( i = 0 ; i < numthreads ; i++ )
	{
		if ( pthread_create(&pthreads[i],NULL, run, (void*)((char*)param+size*i)) != 0 )			
		{
			print_error(__FILE__,(char*)__FUNCTION__,__LINE__,
					    "create %dth thread", i);		
		}
	}
	for ( i = 0 ; i < numthreads ; i++ )
	{
		pthread_join(pthreads[i], NULL);
	}
	memory_free(pthreads);
	pthreads = NULL;
}

void thread_run_with_arraylist(int numthreads, void* (*run)(void*), pArrayList pal )
{
	int i;
	pthread_t* pthreads = (pthread_t*)memory_new( numthreads, sizeof(pthread_t));
	for ( i = 0 ; i < numthreads ; i++ )
	{
		if ( pthread_create(&pthreads[i],NULL, run, pal->array[i] ) != 0 )			
		{
			print_error(__FILE__,(char*)__FUNCTION__,__LINE__,
					    "create %dth thread", i);		
		}
	}
	for ( i = 0 ; i < numthreads ; i++ )
	{
		pthread_join(pthreads[i], NULL);
	}
	memory_free(pthreads);
	pthreads = NULL;
}


void cafe_lambda_set_default(pCafeParam param, double* lambda)
{
	int i;
	pTree ptree = (pTree)param->pcafe;
	for ( i = 0 ; i < ptree->nlist->size; i++ )
	{
		((pCafeNode)ptree->nlist->array[i])->birth_death_probabilities.lambda = lambda[0];
	}
}

void show_sizes(FILE* f, pCafeTree pcafe, family_size_range* range, pCafeFamilyItem pitem, int i)
{
	fprintf(f, ">> %d %d\n", i, pitem->ref);
	fprintf(f, "Root size: %d ~ %d , %d \n",
		pcafe->range.root_min,
		pcafe->range.root_max, pcafe->rfsize);
	fprintf(f, "Family size: %d ~ %d\n", pcafe->range.min, pcafe->range.max);
	fprintf(f, "Root size: %d ~ %d\n", range->root_min, range->root_max);
	fprintf(f, "Family size: %d ~ %d\n", range->min, range->max);
}

void input_values_randomize(input_values *vals, int lambda_len, int mu_len, int k,
	int kfix, double max_branch_length, double *k_weights)
{
	int i,j;
	if (mu_len < 0) {mu_len = 0;}
	if (k > 0) 
	{
		for ( i = 0; i < lambda_len*kfix; i++ )
		{
			vals->parameters[i] = 1.0/max_branch_length * unifrnd();
		}
		int first_mu = lambda_len*kfix;
		for ( i = 0; i < mu_len*kfix; i++ )
		{
			vals->parameters[first_mu+i] = 1.0/max_branch_length * unifrnd();
		}
		double sumofweights = 0;
		for (j = 0; j < k; j++) {
			k_weights[j] = unifrnd();
			sumofweights += k_weights[j];
		}
		for (j = 0; j < k; j++) {
			k_weights[j] = k_weights[j]/sumofweights;
		}
		for (j = 0; j < k-1; j++) {
			vals->parameters[(lambda_len+mu_len)*(kfix)+j] = k_weights[j];
		}
	}
	else {
		for ( i = 0; i < lambda_len; i++ )
		{
			vals->parameters[i] = 1.0/max_branch_length * unifrnd();
		}
		int first_mu = lambda_len;
		for ( i = 0; i < mu_len; i++ )
		{
			vals->parameters[first_mu+i] = 1.0/max_branch_length * unifrnd();
		}
	}
}

double cafe_get_clustered_posterior(pCafeParam param, double *ML, double *MAP, double *prior_rfsize)
{
	int i,j,k;
	double score = 0;
	double** k_likelihoods = NULL;
	double* sumofweights = (double*) memory_new(param->parameterized_k_value, sizeof(double));
	for ( i = 0 ; i < param->pfamily->flist->size ; i++ )
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
		if ( pitem->ref < 0 || pitem->ref == i ) 
		{
			cafe_family_set_size(param->pfamily, pitem, param->pcafe);
			k_likelihoods = cafe_tree_clustered_likelihood(param->pcafe, &cache);		// likelihood of the whole tree = multiplication of likelihood of all nodes
			
			// find the p_z_membership conditioned on the current parameter.
			// it is just proportional to the likelihood of each datapoint in each cluster weighted by the k_weights.
			double sumLikelihood = 0;
			double* MAP_k = memory_new(param->parameterized_k_value, sizeof(double));
			for (k = 0; k < param->parameterized_k_value; k++) {
				
				// get posterior by adding lnPrior to lnLikelihood
				double* posterior = (double*)memory_new(FAMILYSIZEMAX,sizeof(double));
				if(prior_rfsize) {		// prior is a poisson distribution on the root size based on leaves' size
					for(j = 0; j < param->pcafe->rfsize; j++)	// j: root family size
					{
						posterior[j+param->pcafe->range.root_min] = exp(log(k_likelihoods[k][j])+log(prior_rfsize[j]));
					}				
				}
				// Max_on_rootsize k_likelihoods(rootsize)
				MAP_k[k] = __max(posterior,FAMILYSIZEMAX) * param->k_weights[k];			// this part find root size condition with maxlikelihood for each family			
				sumLikelihood += MAP_k[k];
				memory_free(posterior);
				posterior = NULL;
				
			}
			// normalize the ML_k so it becomes a probability
			for (k = 0; k < param->parameterized_k_value; k++) {
				param->p_z_membership[i][k] = MAP_k[k]/sumLikelihood;
				sumofweights[k] += param->p_z_membership[i][k];
			}
			// now since we have the (soft)membership count, we can get the expected logLikelihood given the data and p_z_membership
			// the expected logLikelihood is the weighted sum of loglikelihoods by their soft-membership to each cluster.
			double expectedPosterior = 0;
			for (k = 0; k<param->parameterized_k_value; k++) {
				expectedPosterior += param->p_z_membership[i][k]*(MAP_k[k]);
			}
			MAP[i] = expectedPosterior;
			free(MAP_k);
			

			// what does this do?? I don't know
			if ( pitem->maxlh < 0 )
			{
				int max_k = __maxidx(param->p_z_membership[i],k);
				pitem->maxlh = __maxidx(k_likelihoods[max_k],param->pcafe->rfsize);	
			}
		}
		else
		{
			ML[i] = ML[pitem->ref];
			MAP[i] = MAP[pitem->ref];
			for (k = 0; k < param->parameterized_k_value; k++) {
				param->p_z_membership[i][k] = param->p_z_membership[pitem->ref][k];
				sumofweights[k] += param->p_z_membership[i][k];
			}
			
		}
		if ( MAP[i] == 0 )
		{ 
			show_sizes(stdout, param->pcafe, &param->family_size, pitem, i);
			pString pstr = cafe_tree_string_with_familysize_lambda(param->pcafe);
			fprintf(stderr, "%d: %s\n", i, pstr->buf);
			string_free(pstr);

			score = log(0);
			break;
		}
		score += log(MAP[i]);			// add log-posterior across all families
	}
	for (k = 0; k < param->parameterized_k_value; k++) {
		param->k_weights[k] = sumofweights[k]/param->pfamily->flist->size;
		//fprintf(stdout, "p%d: %f\n", k, param->k_weights[k]);
		/*if (param->k_weights[k] < 2*MIN_DOUBLE) {
		 score = log(0);			// forcing it to be -inf does NOT work, gets stuck in -inf
		 }*/
	}
	memory_free(sumofweights);
	return score;
}

double __lnLGamma(double* palphabeta, void* data)
{
	int i=0;
	double score = 0;
	double alpha = palphabeta[0];
	double beta = palphabeta[1];
	
	pArrayList pal = (pArrayList) data;
	for (i=0; i<pal->size; i++) {
		int* p_x = (int*)pal->array[i];
		int x = *p_x;
		score += log(gampdf((double)x, alpha, beta));
	}
	printf("alpha: %f, beta: %f (Gamma) & Score: %f\n", alpha, beta, score);	
	return -score;
}

double __cafe_cluster_lambda_search(double* parameters, void* args)
{
	int i;
	pCafeParam param = (pCafeParam)args;
	pCafeTree pcafe = (pCafeTree)param->pcafe;
	double score = 0;
	int skip = 0;
	for ( i = 0 ; i < param->num_params; i++ )
	{
		if ( parameters[i] < 0 ) 
		{ 
			skip  = 1;
			score = log(0);
			break;
		}
	}
	if ( !skip )
	{
        cafe_shell_set_lambdas(param, parameters);

		reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, &param->family_size);
		score = cafe_get_clustered_posterior(param, param->ML, param->MAP, param->prior_rfsize);
		cafe_free_birthdeath_cache(pcafe);
		cafe_tree_node_free_clustered_likelihoods(param);
	}
	char buf[STRING_STEP_SIZE];
	buf[0] = '\0';
	string_pchar_join_double(buf,",", param->num_lambdas*(param->parameterized_k_value-param->fixcluster0), parameters );
	fprintf(stdout, "Lambda : %s\n", buf);
	buf[0] = '\0';
	if (param->parameterized_k_value > 0) {
		string_pchar_join_double(buf,",", param->parameterized_k_value, param->k_weights );
		fprintf(stdout, "p : %s\n", buf);
	}
	fprintf(stdout, "Score: %f\n", score);
	fprintf(stdout, ".");
	return -score;
}







extern int chooseln_cache_size;

void reset_birthdeath_cache(pCafeTree tree, int k_value, family_size_range* range)
{
	if (probability_cache)
	{
		birthdeath_cache_array_free(probability_cache);
	}
    probability_cache = cafe_tree_set_birthdeath(tree, MAX(range->max, range->root_max));
}

/**************************************************************************
 * Likelihood ratio test
**************************************************************************/
typedef struct
{
	pCafeParam cafeparam;
	double* maximumPvalues;
	int from;
}LRTParam;

typedef LRTParam* pLRTParam;

pthread_mutex_t mutex_cafe_likelihood = PTHREAD_MUTEX_INITIALIZER;

void* __cafe_likelihood_ratio_test_thread_func(pLRTParam plrt)
{
	int i,b;
	pCafeParam param  = plrt->cafeparam;
	pTree ptree = (pTree)param->pcafe;
	pCafeTree pcafe = cafe_tree_copy(param->pcafe);
	int nnodes = ptree->nlist->size;
	int old_bl;
	int fsize = param->pfamily->flist->size;
#ifdef VERBOSE
	printf("Likelihood ratio test: %d\n", plrt->from );
#endif
	for ( i = plrt->from ; i < fsize ; i += param->num_threads )
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
		if ( pitem->ref >= 0 &&  pitem->ref != i ) continue;
		if (plrt->maximumPvalues[i] > param->pvalue )
		{
			for( b = 0 ; b < nnodes ; b++ ) param->likelihoodRatios[b][i] = -1;
			continue;
		}
		cafe_family_set_size(param->pfamily, pitem, pcafe);
		compute_tree_likelihoods(pcafe);
		double maxlh = __max(get_likelihoods(pcafe), param->pcafe->rfsize );
		for( b = 0 ; b < nnodes ; b++ )
		{
			pPhylogenyNode pnode = (pPhylogenyNode)pcafe->super.nlist->array[b];
			if ( tree_is_root( (pTree)pcafe ,(pTreeNode)pnode) ) 
			{
				param->likelihoodRatios[b][i] = -1;
				continue;
			}
			old_bl = pnode->branchlength;
			struct square_matrix *old_bd = ((pCafeNode)pnode)->birthdeath_matrix;
			double prevlh = -1;
			double nextlh = maxlh;
			while( prevlh < nextlh )
			{
				prevlh = nextlh;
				pnode->branchlength += rint(pnode->branchlength * 0.15);
pthread_mutex_lock( &mutex_cafe_likelihood );
				((pCafeNode)pnode)->birthdeath_matrix = birthdeath_cache_get_matrix(probability_cache, pnode->branchlength, ((pCafeNode)pnode)->birth_death_probabilities.lambda,  ((pCafeNode)pnode)->birth_death_probabilities.mu );
pthread_mutex_unlock( &mutex_cafe_likelihood);
				compute_tree_likelihoods(pcafe);
				nextlh = __max(get_likelihoods(pcafe), param->pcafe->rfsize );
			}
			param->likelihoodRatios[b][i] = (prevlh == maxlh) ? 1 : 1 - chi2cdf( 2*(log(prevlh) - log(maxlh)), 1);
		//	param->likelihoodRatios[b][i] = (prevlh == maxlh) ? 1 : prevlh / maxlh ;
			pnode->branchlength = old_bl;
			((pCafeNode)pnode)->birthdeath_matrix = old_bd;
		}
	}
	cafe_tree_free(pcafe);
	return (NULL);
}

void cafe_likelihood_ratio_test(pCafeParam param, double *maximumPvalues)
{
	cafe_log(param,"Running Likelihood Ratio Test....\n");

	pTree ptree = (pTree)param->pcafe;
	int nnodes = ptree->nlist->size;

	pLRTParam ptparam = (pLRTParam)memory_new(param->num_threads,sizeof(LRTParam));
	int nrows = param->pfamily->flist->size;

	param->likelihoodRatios = (double**)memory_new_2dim(nnodes,nrows,sizeof(double));

	int i, b;
	for ( i = 0 ; i < param->num_threads; i++)
	{
		ptparam[i].cafeparam = param;
		ptparam[i].from = i;
		ptparam[i].maximumPvalues = maximumPvalues;
        __cafe_likelihood_ratio_test_thread_func(ptparam);
    }
	for( i = 0 ; i < nrows ; i++ )
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
		if ( pitem->ref < 0 || pitem->ref == i ) continue;
		for( b = 0 ; b < nnodes ; b++ )
		{
			param->likelihoodRatios[b][i] = param->likelihoodRatios[b][pitem->ref];
		}
		continue;
	}
	memory_free(ptparam);	
	ptparam = NULL;				
	cafe_log( param , "Done : Likelihood Ratio test\n" );
}
