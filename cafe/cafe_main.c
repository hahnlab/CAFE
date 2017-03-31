/*! \mainpage Index Page
*
* \section intro_sec Introduction
*
* gene family evolution, applied in the software package, CAFE. 
Application of this method to data from multiple whole genomes 
of many groups is revealing remarkable patterns of gene gain 
and loss. Other approaches to studying this question have 
involved the analysis of gene movement among chromosomes 
(especially sex chromosomes), the discovery of polymorphic 
copy-number variants under local selection, and even new 
methods for carrying out genome assembly to more accurately 
estimate gene numbers.
*
* \section install_sec Installation
*
* \subsection step1 Download and make
*
* \subsection Program structure
*
* The list of available commands are stored in the #cafe_cmd
variable in cafe_shell.c .
* 
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
	tree->rootfamilysizes[0] = range->root_min;
	tree->rootfamilysizes[1] = range->root_max;
	tree->familysizes[0] = range->min;
	tree->familysizes[1] = range->max;
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
		pcafe->rootfamilysizes[0],
		pcafe->rootfamilysizes[1], pcafe->rfsize);
	fprintf(f, "Family size: %d ~ %d\n", pcafe->familysizes[0], pcafe->familysizes[1]);
	fprintf(f, "Root size: %d ~ %d\n", range->root_min, range->root_max);
	fprintf(f, "Family size: %d ~ %d\n", range->min, range->max);
}

double cafe_get_posterior(pCafeFamily pfamily, pCafeTree pcafe, family_size_range*range, double *ML, double *MAP, double *prior_rfsize, int quiet)
{
	int i, j;
	double score = 0;
	double* likelihood = NULL;
	for ( i = 0 ; i < pfamily->flist->size ; i++ )	// i: family index
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)pfamily->flist->array[i];
		if ( pitem->ref < 0 || pitem->ref == i ) 
		{
			cafe_family_set_size(pfamily, i, pcafe);	// this part is just setting the leave counts.
			compute_tree_likelihoods(pcafe);
			likelihood = get_likelihoods(pcafe);		// likelihood of the whole tree = multiplication of likelihood of all nodes
			ML[i] = __max(likelihood, pcafe->rfsize);			// this part find root size condition with maxlikelihood for each family			
			if ( pitem->maxlh < 0 )
			{
				pitem->maxlh = __maxidx(likelihood, pcafe->rfsize);	
			}
			// get posterior by adding lnPrior to lnLikelihood
			double* posterior = (double*)memory_new(pcafe->size_of_factor,sizeof(double));
			if(prior_rfsize) {		// prior is a poisson distribution on the root size based on leaves' size
				for(j = 0; j < pcafe->rfsize; j++)	// j: root family size
				{
					// likelihood and posterior both starts from 1 instead of 0 
					posterior[j] = exp(log(likelihood[j])+log(prior_rfsize[j]));	//prior_rfsize also starts from 1
				}				
			}
            else {
                fprintf(stderr,"ERROR: empirical posterior not defined.\n");      
                return -1;
            }
			MAP[i] = __max(posterior, pcafe->rfsize);			// this part find root size condition with maxlikelihood for each family			
			memory_free(posterior);
			posterior = NULL;
		}
		else
		{
			ML[i] = ML[pitem->ref];
			MAP[i] = MAP[pitem->ref];
		}
		if ( ML[i] == 0 )
		{ 
			if (!quiet)
			{ 
				show_sizes(stdout, pcafe, range, pitem, i);
				pString pstr = cafe_tree_string_with_familysize_lambda(pcafe);
				fprintf(stderr, "%d: %s\n", i, pstr->buf );
				string_free(pstr);
			}
			score = log(0);
			break;
		}
		score += log(MAP[i]);			// add log-posterior across all families
	}
	return score;
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

double cafe_get_clustered_posterior(pCafeParam param)
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
			cafe_family_set_size(param->pfamily, i, param->pcafe);
			k_likelihoods = cafe_tree_clustered_likelihood(param->pcafe);		// likelihood of the whole tree = multiplication of likelihood of all nodes
			
			// find the p_z_membership conditioned on the current parameter.
			// it is just proportional to the likelihood of each datapoint in each cluster weighted by the k_weights.
			double sumLikelihood = 0;
			double* MAP_k = memory_new(param->parameterized_k_value, sizeof(double));
			for (k = 0; k < param->parameterized_k_value; k++) {
				
				// get posterior by adding lnPrior to lnLikelihood
				double* posterior = (double*)memory_new(FAMILYSIZEMAX,sizeof(double));
				if(param->prior_rfsize) {		// prior is a poisson distribution on the root size based on leaves' size
					for(j = 0; j < param->pcafe->rfsize; j++)	// j: root family size
					{
						posterior[j+param->pcafe->rootfamilysizes[0]] = exp(log(k_likelihoods[k][j])+log(param->prior_rfsize[j]));
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
			param->MAP[i] = expectedPosterior;
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
			param->ML[i] = param->ML[pitem->ref];
			param->MAP[i] = param->MAP[pitem->ref];
			for (k = 0; k < param->parameterized_k_value; k++) {
				param->p_z_membership[i][k] = param->p_z_membership[pitem->ref][k];
				sumofweights[k] += param->p_z_membership[i][k];
			}
			
		}
		if ( param->MAP[i] == 0 )
		{ 
			show_sizes(stdout, param->pcafe, &param->family_size, pitem, i);
			pString pstr = cafe_tree_string_with_familysize_lambda(param->pcafe);
			fprintf(stderr, "%d: %s\n", i, pstr->buf);
			string_free(pstr);

			score = log(0);
			break;
		}
		score += log(param->MAP[i]);			// add log-posterior across all families
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

double __lnLPoisson(double* plambda, void* data)
{
	int i=0;
	double score = 0;
	double lambda = plambda[0];
	pArrayList pal = (pArrayList) data;
	for (i=0; i<pal->size; i++) {
		int* p_x = (int*)pal->array[i];
		int x = *p_x;
		double ll = poisspdf((double)x, lambda);
		if (isnan(ll)) {
			ll = 0;
		}
		score += log(ll);
	}
	//printf("lambda: %f (Poisson) & Score: %f\n", lambda, score);	
	return -score;
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




double cafe_set_prior_rfsize_poisson_lambda(pCafeParam param, double* lambda)
{
	int i;
	// calculate the prior probability for a range of root sizes.
	if ( param->prior_rfsize )
	{
		memory_free(param->prior_rfsize);
		param->prior_rfsize = NULL;
	}
	param->prior_rfsize = (double *) memory_new(FAMILYSIZEMAX, sizeof(double));
	for( i=0; i<FAMILYSIZEMAX; i++) {
		//param->prior_rfsize[i] = poisspdf(param->pcafe->rootfamilysizes[0]+i, parameters[0]);					// poisson
		param->prior_rfsize[i] = poisspdf(param->pcafe->rootfamilysizes[0]-1+i, lambda[0]);					// shifted poisson
		//param->prior_rfsize[i] = gampdf(param->pcafe->rootfamilysizes[0]+i, parameters[0], parameters[1]);	// gamma
	}	
	return 0;
}

/// set empirical prior on rootsize based on the assumption that rootsize follows leaf size distribution
double cafe_set_prior_rfsize_empirical(pCafeParam param)
{
	int i=0;
	int idx = 0;

	// estimate the distribution of family size based on observed leaf counts.
	// first collect all leaves sizes into an ArrayList.
	pArrayList pLeavesSize = arraylist_new(param->pfamily->flist->size*param->pfamily->num_species);	
	for( idx=0; idx<param->pfamily->flist->size; idx++) {
		pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[idx];
		for ( i = 0 ; i < param->pfamily->num_species ; i++ )
		{
			if( param->pfamily->index[i] < 0 ) continue;
			int* leafcnt = memory_new(1, sizeof(int));
			memcpy(leafcnt, &(pitem->count[i]), sizeof(int));
			if (*leafcnt > 0) {		// ignore the zero counts ( we condition that rootsize is at least one )
				*leafcnt = (*leafcnt) - 1;
				arraylist_add(pLeavesSize, (void*)leafcnt);
			}
		}
	}

	// now estimate parameter based on data and distribution (poisson or gamma). 
	pFMinSearch pfm;
	int num_params = 1;
	pfm = fminsearch_new_with_eq(__lnLPoisson,num_params,pLeavesSize);
	//int num_params = 2;
	//pfm = fminsearch_new_with_eq(__lnLGamma,num_params,pLeavesSize);
	pfm->tolx = 1e-6;
	pfm->tolf = 1e-6;
	double* parameters = memory_new(num_params, sizeof(double));
	for ( i = 0; i < num_params; i++ ) parameters[i] = unifrnd();
	fminsearch_min(pfm, parameters);
	double *re = fminsearch_get_minX(pfm);
	for ( i = 0; i < num_params; i++ ) parameters[i] = re[i];
	cafe_log(param,"Empirical Prior Estimation Result: %d\n", pfm->iters );
	cafe_log(param,"Poisson lambda: %f & Score: %f\n", parameters[0], *pfm->fv);	
	double *prior_poisson_lambda = memory_new_with_init(num_params, sizeof(double), (void*) parameters);
	//cafe_log(param,"Gamma alpha: %f, beta: %f & Score: %f\n", parameters[0], parameters[1], *pfm->fv);	
	
	// set rfsize based on estimated prior
	cafe_set_prior_rfsize_poisson_lambda(param, prior_poisson_lambda);

	// clean
	fminsearch_free(pfm);
	arraylist_free(pLeavesSize, NULL);
	memory_free(parameters);
	return 0;
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
		param->param_set_func(param,parameters);

		reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, &param->family_size);
		score = cafe_get_clustered_posterior(param);
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




double __cafe_best_lambda_search(double* plambda, void* args)
{
	int i;
	pCafeParam param = (pCafeParam)args;
	pCafeTree pcafe = (pCafeTree)param->pcafe;
	double score = 0;
	int skip = 0;
	for ( i = 0 ; i < param->num_lambdas ; i++ )
	{
		if ( plambda[i] < 0 ) 
		{ 
			skip  = 1;
			score = log(0);
			break;
		}
	}
	if ( !skip )
	{
		param->param_set_func(param,plambda);

		reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, &param->family_size);
        score = cafe_get_posterior(param->pfamily, param->pcafe, &param->family_size, param->ML, param->MAP, param->prior_rfsize, param->quiet);
		cafe_free_birthdeath_cache(pcafe);
	}
	char buf[STRING_STEP_SIZE];
	buf[0] = '\0';
	string_pchar_join_double(buf,",", param->num_lambdas, plambda );
	cafe_log(param,"Lambda : %s & Score: %f\n", buf, score);
	cafe_log(param, ".");
	return -score;
}


extern int chooseln_cache_size;

double* cafe_best_lambda_by_fminsearch(pCafeParam param, int lambda_len, int k )
{
	int i,j;
	int max_runs = 10;
	double* scores = memory_new(max_runs, sizeof(double));
	int converged = 0;
	int runs = 0;
	
	do
	{
		
		if ( param->num_params > 0 )
		{
			int kfix = k - param->fixcluster0;
			double max_branch_length = param->max_branch_length;
			double *k_weights = param->k_weights;
			input_values_randomize(&param->input, param->num_lambdas, param->num_mus, param->parameterized_k_value,
				kfix, max_branch_length, k_weights);
		}
		
		copy_range_to_tree(param->pcafe, &param->family_size);
		
		pFMinSearch pfm;
		if (k > 0) {
			pfm = fminsearch_new_with_eq(__cafe_cluster_lambda_search, param->num_params, param);
			pfm->tolx = 1e-5;
			pfm->tolf = 1e-5;
		}
		else {
			pfm = fminsearch_new_with_eq(__cafe_best_lambda_search,lambda_len,param);
			pfm->tolx = 1e-6;
			pfm->tolf = 1e-6;
		}
		fminsearch_min(pfm, param->input.parameters);
		double *re = fminsearch_get_minX(pfm);
		for ( i = 0 ; i < param->num_params ; i++ ) param->input.parameters[i] = re[i];
        
		double current_p = param->input.parameters[(lambda_len)*(k-param->fixcluster0)];
		double prev_p;
		if (k>0) 
		{
			do {
				double* sumofweights = (double*) memory_new(param->parameterized_k_value, sizeof(double));
				for ( i = 0 ; i < param->pfamily->flist->size ; i++ ) {
					for (j = 0; j<k; j++) {
						sumofweights[j] += param->p_z_membership[i][j];
					}
				}
				for (j = 0; j<k-1; j++) {
					param->input.parameters[(lambda_len)*(k-param->fixcluster0)+j] = sumofweights[j]/param->pfamily->flist->size;
				}
				memory_free(sumofweights);

				fminsearch_min(pfm, param->input.parameters);
				
				double *re = fminsearch_get_minX(pfm);
				for ( i = 0 ; i < param->num_params ; i++ ) param->input.parameters[i] = re[i];
				
				prev_p = current_p;
				current_p = param->input.parameters[(lambda_len)*(k-param->fixcluster0)];
			} while (current_p - prev_p > pfm->tolx);
		}
		
		cafe_log(param, "\n");
		cafe_log(param,"Lambda Search Result: %d\n", pfm->iters );
		if (k > 0) {
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			if (param->fixcluster0) {
				strncat(buf, "0,", 2);
				string_pchar_join_double(buf,",", param->num_lambdas*(param->parameterized_k_value-param->fixcluster0), param->input.parameters );
			}
			else {
				string_pchar_join_double(buf,",", param->num_lambdas*param->parameterized_k_value, param->input.parameters );
			}
			cafe_log(param,"Lambda : %s\n", buf);
			buf[0] = '\0';
			if (param->parameterized_k_value > 0) {
				string_pchar_join_double(buf,",", param->parameterized_k_value, param->k_weights );
				cafe_log(param, "p : %s\n", buf);
				cafe_log(param, "p0 : %f\n", param->input.parameters[param->num_lambdas*(param->parameterized_k_value-param->fixcluster0)+0]);
			}
			cafe_log(param, "Score: %f\n", *pfm->fv);
		}
		else {
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			string_pchar_join_double(buf,",", param->num_lambdas, param->input.parameters );
			cafe_log(param,"Lambda : %s & Score: %f\n", buf, *pfm->fv);
		}
		if (runs > 0) {
			double minscore = __min(scores, runs);
			if (abs(minscore - (*pfm->fv)) < 10*pfm->tolf) {
				converged = 1;
			}
		}
		scores[runs] = *pfm->fv;
		fminsearch_free(pfm);
		
		copy_range_to_tree(param->pcafe, &param->family_size);
		
		runs++;

	} while (param->checkconv && !converged && runs<max_runs); 
		
//	string_free(pstr);
	if (param->checkconv) {
		if (converged) {
			cafe_log(param,"score converged in %d runs.\n", runs);
		}
		else {
			cafe_log(param,"score failed to converge in %d runs.\n", max_runs);
		}
	}
	memory_free(scores);
	return param->input.parameters;
}

void reset_birthdeath_cache(pCafeTree tree, int k_value, family_size_range* range)
{
	if (probability_cache)
	{
		birthdeath_cache_array_free(probability_cache);
	}
	probability_cache = birthdeath_cache_init(MAX(range->max, range->root_max));
	cafe_tree_set_birthdeath(tree, probability_cache);
}

double __cafe_each_best_lambda_search(double* plambda, void* args)
{
	int i;
	pCafeParam param = (pCafeParam)args;
	pCafeTree pcafe = (pCafeTree)param->pcafe;
	double score = 0;
	int skip = 0;
	for ( i = 0 ; i < param->num_lambdas ; i++ )
	{
		if ( plambda[i] < 0 ) 
		{ 
			skip  = 1;
			score = log(0);
			break;
		}
	}

	if ( !skip )
	{
		param->param_set_func(param,plambda);

		reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, &param->family_size);
		compute_tree_likelihoods(pcafe);
		double* likelihood = get_likelihoods(pcafe);
		score = log(__max(likelihood,pcafe->rfsize));
		cafe_free_birthdeath_cache(pcafe);
		probability_cache = NULL;
	}

	char buf[STRING_STEP_SIZE];
	buf[0] = '\0';
	string_pchar_join_double(buf,",", param->num_lambdas, plambda );
	cafe_log(param,"\tLambda : %s & Score: %f\n", buf, score);
	cafe_log(param, "\n");
	return -score;

}

double* cafe_each_best_lambda_by_fminsearch(pCafeParam param, int lambda_len )
{
	double* old_lambda = NULL;
	if( param->lambda ) 
	{
		old_lambda = param->lambda;
	}
	param->lambda = (double*)memory_new(lambda_len, sizeof(double));
	param->num_lambdas = lambda_len;

	family_size_range temp_range = param->family_size;

	int i, j;
	for ( i = 0 ; i < lambda_len ; i++ )
	{
		param->lambda[i] = 0.5/param->max_branch_length;
	}
	pFMinSearch pfm = fminsearch_new_with_eq(__cafe_each_best_lambda_search,lambda_len,param);
	pfm->tolx = 1e-6;
	pfm->tolf = 1e-6;
	int fsize = param->pfamily->flist->size;
	for ( i = 0 ; i < param->pfamily->flist->size ; i++ )
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
		if ( pitem->ref >= 0 && pitem->ref != i )
		{
			pCafeFamilyItem pref = (pCafeFamilyItem)param->pfamily->flist->array[pitem->ref];
			pitem->lambda = pref->lambda;
			pitem->mu = pref->mu;
			param->param_set_func(param,pitem->lambda);

			cafe_log(param,"%s: Lambda Search Result of %d/%d in %d iteration \n", pitem->id, i+1, fsize, pfm->iters );
			pString pstr = cafe_tree_string_with_familysize_lambda(param->pcafe);
			cafe_log(param,"%s: %s\n", pitem->id, pstr->buf );
			string_free(pstr);
			continue;
		}

		cafe_family_set_size_with_family_forced(param->pfamily,i,param->pcafe);

		param->family_size.root_min = param->pcafe->rootfamilysizes[0];
		param->family_size.root_max = param->pcafe->rootfamilysizes[1];
		param->family_size.min = param->pcafe->familysizes[0];
		param->family_size.max = param->pcafe->familysizes[1];

		cafe_log(param,"%s:\n", pitem->id );
		
		fminsearch_min(pfm, param->lambda );

/*
		printf("%d %d:", param->rootfamily_sizes[1], param->family_sizes[1] );
		for ( j = 0 ; j < param->pcafe->super.nlist->size ; j+=2 )
		{
			pCafeNode pnode = (pCafeNode)param->pcafe->super.nlist->array[j];
			printf(" %d",  pnode->familysize );
		}
		printf("\n");
*/

		double *re = fminsearch_get_minX(pfm);
		if ( pitem->lambda ) {memory_free( pitem->lambda ); pitem->lambda = NULL; }
		if ( pitem->mu ) {memory_free( pitem->mu ); pitem->mu = NULL;}
		pitem->lambda = (double*) memory_new(lambda_len, sizeof(double));
		pitem->mu = (double*) memory_new(lambda_len, sizeof(double));
		int lambda_check = 0;
		for ( j = 0 ; j < lambda_len ; j++ ) 
		{
			pitem->lambda[j] = re[j];
			//pitem->mu[j] = re[j];
			double a = re[j] * param->max_branch_length;
//			printf("%f\n", a);
			if ( a >= 0.5 || fabs(a-0.5) < 1e-3 )
			{
				lambda_check = 1;	
			}	
		}
		param->param_set_func(param,re);

		cafe_log(param,"Lambda Search Result of %d/%d in %d iteration \n", i+1, fsize, pfm->iters );
		if ( lambda_check )
		{
			cafe_log(param,"Caution : at least one lambda near boundary\n" );
		}
		pString pstr = cafe_tree_string_with_familysize_lambda(param->pcafe);
		if ( lambda_check )
		{
			cafe_log(param,"@@ ");
		}
		cafe_log(param,"%s\n", pstr->buf );
		string_free(pstr);
	}
	fminsearch_free(pfm);

	copy_range_to_tree(param->pcafe, &temp_range);

	param->family_size = temp_range;

	memory_free(param->lambda);
	param->lambda = old_lambda;
	return param->lambda;
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

void* __cafe_likelihood_ratio_test_thread_func(void* ptr)
{
	int i,b;
	pLRTParam plrt = (pLRTParam)ptr;
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
		cafe_family_set_size(param->pfamily,i, pcafe);
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
	}
	thread_run(param->num_threads, __cafe_likelihood_ratio_test_thread_func, ptparam, sizeof(LRTParam));
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

/*******************************************************************************
 *	Cafe Parameter
 *******************************************************************************/

pCafeParam cafe_copy_parameters(pCafeParam psrc)
{
	pCafeParam param = (pCafeParam)memory_new(1, sizeof(CafeParam));
	memcpy( param, psrc, sizeof(CafeParam));
	param->lambda = NULL;
	param->num_lambdas = 0;
	param->pcafe = cafe_tree_copy(psrc->pcafe);

//	param->viterbi.viterbiPvalues = NULL;
//	param->viterbi.expandRemainDecrease = NULL;
//	param->viterbi.viterbiNodeFamilysizes = NULL;
//	param->viterbi.maximumPvalues = NULL;
//	param->viterbi.averageExpansion = NULL;
//	param->viterbi.cutPvalues = NULL;

	return param;
}

void cafe_free_copy_parameters(pCafeParam param)
{
	//memory_free(param->branchlengths_sorted);	
	//param->branchlengths_sorted = NULL;
	cafe_tree_free(param->pcafe);
	memory_free(param);
	param = NULL;
}
