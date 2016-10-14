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

void cafe_log(pCafeParam param, const char* msg, ... )
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

void cafe_add_birthdeath_cache(pCafeParam param, pCafeTree pcafe )
{
	int i;
	for ( i = 0 ; i < pcafe->super.nlist->size ; i++ )
	{
		pPhylogenyNode pnode = (pPhylogenyNode)pcafe->super.nlist->array[i];
		if ( pnode->branchlength > 0 )
		{
			birthdeath_cache_get_matrix(pcafe->pbdc_array, pnode->branchlength , ((pCafeNode)pnode)->lambda, ((pCafeNode)pnode)->mu );
		}
	}
	cafe_tree_set_birthdeath(param->pcafe);
}


void cafe_free_birthdeath_cache(pCafeTree pcafe)
{
	birthdeath_cache_array_free(pcafe->pbdc_array);
	pcafe->pbdc_array = NULL;
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
		((pCafeNode)ptree->nlist->array[i])->lambda = lambda[0];
	}
}

/**************************************************************************
 * Lambda 
**************************************************************************/
double __cafe_best_lambda_search(double* plambda, void* args);

pGMatrix cafe_lambda_distribution(pCafeParam param, int numrange, double** range )
{
	int i, j;
	int* size = (int*)memory_new(numrange,sizeof(int));
	double* plambda = (double*)memory_new(numrange,sizeof(double)); 
	int* idx = (int*)memory_new(numrange, sizeof(int));
	for ( i = 0 ; i < numrange; i++ )
	{
		size[i] = 1 + rint((range[i][2] - range[i][0])/range[i][1]);
	}
	pGMatrix pgm = gmatrix_double_new(numrange,size);

	for ( i = 0 ; i < pgm->num_elements; i++ )
	{
		gmatrix_dim_index(pgm,i,idx);
		for ( j = 0 ; j < numrange ; j++ )
		{
			plambda[j] = range[j][1] * idx[j] + range[j][0];
		}
		double v = -__cafe_best_lambda_search(plambda,(void*)param);
		gmatrix_double_set_with_index(pgm, v, i );
		if ( -v > 1e300 )
		{
			for ( j = 0 ; j < param->pfamily->flist->size ; j++ )
			{
				pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[j];
				pitem->maxlh = -1;
			}
		}
	}

	memory_free(idx);
	idx = NULL;
	memory_free(size);
	size = NULL;
	memory_free(plambda);
	plambda = NULL;
	return pgm;
}



void show_sizes(FILE* f, pCafeParam param, pCafeFamilyItem pitem, int i)
{
	fprintf(f, ">> %d %d\n", i, pitem->ref);
	fprintf(f, "Root size: %d ~ %d , %d \n",
		param->pcafe->rootfamilysizes[0],
		param->pcafe->rootfamilysizes[1], param->pcafe->rfsize);
	fprintf(f, "Family size: %d ~ %d\n", param->pcafe->familysizes[0], param->pcafe->familysizes[1]);
	fprintf(f, "Root size: %d ~ %d\n", param->rootfamily_sizes[0], param->rootfamily_sizes[1]);
	fprintf(f, "Family size: %d ~ %d\n", param->family_sizes[0], param->family_sizes[1]);
}

double cafe_get_posterior(pCafeParam param)
{
	int i, j;
	double score = 0;
	double* likelihood = NULL;
	for ( i = 0 ; i < param->pfamily->flist->size ; i++ )	// i: family index
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
		if ( pitem->ref < 0 || pitem->ref == i ) 
		{
			cafe_family_set_size(param->pfamily, i, param->pcafe);	// this part is just setting the leave counts.
			likelihood = cafe_tree_likelihood(param->pcafe);		// likelihood of the whole tree = multiplication of likelihood of all nodes
			param->ML[i] = __max(likelihood,param->pcafe->rfsize);			// this part find root size condition with maxlikelihood for each family			
			if ( pitem->maxlh < 0 )
			{
				pitem->maxlh = __maxidx(likelihood,param->pcafe->rfsize);	
			}
			// get posterior by adding lnPrior to lnLikelihood
			double* posterior = (double*)memory_new(param->pcafe->size_of_factor,sizeof(double));
			if (param->prior_rfsize_by_family) {		// prior is set by birth=death
				for(j = 0; j < param->pcafe->rfsize; j++)	// j: root family size
				{
					posterior[j] = exp(log(likelihood[j])+log(param->prior_rfsize_by_family[i][j]));
				}
			}
			else if(param->prior_rfsize) {		// prior is a poisson distribution on the root size based on leaves' size
				for(j = 0; j < param->pcafe->rfsize; j++)	// j: root family size
				{
					// likelihood and posterior both starts from 1 instead of 0 
					posterior[j] = exp(log(likelihood[j])+log(param->prior_rfsize[j]));	//prior_rfsize also starts from 1
				}				
			}
            else {
                fprintf(stderr,"ERROR: empirical posterior not defined.\n");      
                return -1;
            }
			param->MAP[i] = __max(posterior,param->pcafe->rfsize);			// this part find root size condition with maxlikelihood for each family			
			memory_free(posterior);
			posterior = NULL;
		}
		else
		{
			param->ML[i] = param->ML[pitem->ref];
			param->MAP[i] = param->MAP[pitem->ref];
		}
		if ( param->ML[i] == 0 )
		{ 
			if (!param->quiet)
			{ 
				show_sizes(stdout, param, pitem, i);
				pString pstr = cafe_tree_string_with_familysize_lambda(param->pcafe);
				fprintf(stderr, "%d: %s\n", i, pstr->buf );
				string_free(pstr);
			}
			score = log(0);
			break;
		}
		score += log(param->MAP[i]);			// add log-posterior across all families
	}
	return score;
}

void __cafe_randomize_cluster_parameters(pCafeParam param, int lambda_len, int mu_len, int k) 
{
	int i,j;
	if (mu_len < 0) {mu_len = 0;}
	if (k > 0) 
	{
		for ( i = 0; i < lambda_len*(k-param->fixcluster0); i++ )
		{
			param->parameters[i] = 1.0/param->max_branch_length * unifrnd();
		}
		for ( i = 0; i < mu_len*(k-param->fixcluster0); i++ )
		{
            int idx = lambda_len*(k-param->fixcluster0)+i;
			param->parameters[idx] = 1.0/param->max_branch_length * unifrnd();
		}
		double sumofweights = 0;
		for (j = 0; j < k; j++) {
			param->k_weights[j] = unifrnd();
			sumofweights += param->k_weights[j];
		}
		for (j = 0; j < k; j++) {
			param->k_weights[j] = param->k_weights[j]/sumofweights;
		}
		for (j = 0; j < k-1; j++) {
			param->parameters[(lambda_len+mu_len)*(k-param->fixcluster0)+j] = param->k_weights[j];
		}
	}
	else {
		for ( i = 0; i < lambda_len; i++ )
		{
			param->parameters[i] = 1.0/param->max_branch_length * unifrnd();
		}
		for ( i = 0; i < mu_len; i++ )
		{
            int idx = lambda_len+i;
			param->parameters[idx] = 1.0/param->max_branch_length * unifrnd();
		}
	}
}

void __cafe_scaleup_cluster_parameters(pCafeParam param, int lambda_len, int mu_len, int k) 
{
	int i;
	if (mu_len < 0) {mu_len = 0;}
	if (k > 0) 
	{
		for ( i = 0; i < lambda_len*(k-param->fixcluster0); i++ )
		{
            param->parameters[i] = param->parameters[i]*param->sum_branch_length; //scale by sum_branch_length;
		}
		for ( i = 0; i < mu_len*(k-param->fixcluster0); i++ )
		{
            int idx = lambda_len*(k-param->fixcluster0)+i;
            param->parameters[idx] = param->parameters[idx]*param->sum_branch_length; //scale by sum_branch_length;
		}
	}
	else {
		for ( i = 0; i < lambda_len; i++ )
		{
            param->parameters[i] = param->parameters[i]*param->sum_branch_length; //scale by sum_branch_length;
		}
		for ( i = 0; i < mu_len; i++ )
		{
            int idx = lambda_len+i;
            param->parameters[idx] = param->parameters[idx]*param->sum_branch_length; //scale by sum_branch_length;
		}
	}
}

void __cafe_scaledown_cluster_parameters(pCafeParam param, int lambda_len, int mu_len, int k) 
{
	int i;
	if (mu_len < 0) {mu_len = 0;}
	if (k > 0) 
	{
		for ( i = 0; i < lambda_len*(k-param->fixcluster0); i++ )
		{
            param->parameters[i] = param->parameters[i]/param->sum_branch_length; //scale by sum_branch_length;
		}
		for ( i = 0; i < mu_len*(k-param->fixcluster0); i++ )
		{
            int idx = lambda_len*(k-param->fixcluster0)+i;
            param->parameters[idx] = param->parameters[idx]/param->sum_branch_length; //scale by sum_branch_length;
		}
	}
	else {
		for ( i = 0; i < lambda_len; i++ )
		{
            param->parameters[i] = param->parameters[i]/param->sum_branch_length; //scale by sum_branch_length;
		}
		for ( i = 0; i < mu_len; i++ )
		{
            int idx = lambda_len+i;
            param->parameters[idx] = param->parameters[idx]/param->sum_branch_length; //scale by sum_branch_length;
		}
	}
}


double cafe_get_clustered_posterior(pCafeParam param)
{
	int i,j,k;
	double score = 0;
	double** k_likelihoods = NULL;
	double* sumofweights = (double*) memory_new(param->k, sizeof(double));
	
	
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
			double* MAP_k = memory_new(param->k, sizeof(double));
			for (k = 0; k < param->k; k++) {
				
				// get posterior by adding lnPrior to lnLikelihood
				double* posterior = (double*)memory_new(FAMILYSIZEMAX,sizeof(double));
				if (param->prior_rfsize_by_family) {		// prior is set by birth=death
					for(j = 0; j < param->pcafe->rfsize; j++)	// j: root family size
					{
						posterior[j+param->pcafe->rootfamilysizes[0]] = exp(log(k_likelihoods[k][j])+log(param->prior_rfsize_by_family[i][j]));
					}
				}
				else if(param->prior_rfsize) {		// prior is a poisson distribution on the root size based on leaves' size
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
			for (k = 0; k < param->k; k++) {
				param->p_z_membership[i][k] = MAP_k[k]/sumLikelihood;
				sumofweights[k] += param->p_z_membership[i][k];
			}
			// now since we have the (soft)membership count, we can get the expected logLikelihood given the data and p_z_membership
			// the expected logLikelihood is the weighted sum of loglikelihoods by their soft-membership to each cluster.
			double expectedPosterior = 0;
			for (k = 0; k<param->k; k++) {
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
			for (k = 0; k < param->k; k++) {
				param->p_z_membership[i][k] = param->p_z_membership[pitem->ref][k];
				sumofweights[k] += param->p_z_membership[i][k];
			}
			
		}
		if ( param->MAP[i] == 0 )
		{ 
			show_sizes(stdout, param, pitem, i);
			pString pstr = cafe_tree_string_with_familysize_lambda(param->pcafe);
			fprintf(stderr, "%d: %s\n", i, pstr->buf);
			string_free(pstr);

			score = log(0);
			break;
		}
		score += log(param->MAP[i]);			// add log-posterior across all families
	}
	for (k = 0; k < param->k; k++) {
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
	param->prior_poisson_lambda = memory_new_with_init(num_params, sizeof(double), (void*) parameters);
	//cafe_log(param,"Gamma alpha: %f, beta: %f & Score: %f\n", parameters[0], parameters[1], *pfm->fv);	
	
	// set rfsize based on estimated prior
	cafe_set_prior_rfsize_poisson_lambda(param, param->prior_poisson_lambda);

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
		cafe_set_birthdeath_cache_thread(param);
//		if (param->posterior) {
			score = cafe_get_clustered_posterior(param);
//		}
//		else {
//			score = cafe_get_clustered_likelihood(param);
//		}
		cafe_free_birthdeath_cache(pcafe);
		cafe_tree_node_free_clustered_likelihoods(param);
	}
	char buf[STRING_STEP_SIZE];
	buf[0] = '\0';
	string_pchar_join_double(buf,",", param->num_lambdas*(param->k-param->fixcluster0), parameters );
	fprintf(stdout, "Lambda : %s\n", buf);
	buf[0] = '\0';
	if (param->k > 0) {
		string_pchar_join_double(buf,",", param->k, param->k_weights );
		fprintf(stdout, "p : %s\n", buf);
	}
	fprintf(stdout, "Score: %f\n", score);
	fprintf(stdout, ".");
	return -score;
}




double __cafe_cluster_lambda_mu_search(double* parameters, void* args)
{
	int i;
	pCafeParam param = (pCafeParam)args;
	pCafeTree pcafe = (pCafeTree)param->pcafe;
	double score = 0;
	int skip = 0;
	for ( i = 0 ; i < param->num_params ; i++ )
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
		cafe_set_birthdeath_cache_thread(param);
//		if (param->posterior) {
			score = cafe_get_clustered_posterior(param);
//		}
//		else {
//			score = cafe_get_clustered_likelihood(param);
//		}
		cafe_free_birthdeath_cache(pcafe);
		cafe_tree_node_free_clustered_likelihoods(param);
	}
	char buf[STRING_STEP_SIZE];
	buf[0] = '\0';
	for( i=0; i<param->num_lambdas; i++) {
	string_pchar_join_double(buf,",", (param->k-param->fixcluster0), &parameters[i*(param->k-param->fixcluster0)] );
	fprintf(stdout, "Lambda branch %d: %s\n", i, buf);
	buf[0] = '\0';
	}
	for (i=0; i<param->num_mus; i++) {
	string_pchar_join_double(buf,",", (param->k-param->fixcluster0), &parameters[param->num_lambdas*(param->k-param->fixcluster0)+i*(param->k-param->fixcluster0)]);
	fprintf(stdout, "Mu branch %d: %s \n", i, buf);
	buf[0] = '\0';
	}
	if (param->k > 0) {
		string_pchar_join_double(buf,",", param->k, param->k_weights );
		fprintf(stdout, "p : %s\n", buf);
	}
	//cafe_log(param, "Score: %f\n", score);
	fprintf(stdout, ".");
	return -score;
}



// need to define a function with lambda and mu. 
// this function is provided as the equation to fmin search.
// also need to make a new param_set_func that includes mu.

double __cafe_best_lambda_mu_search(double* parameters, void* args)
{
	int i;
	pCafeParam param = (pCafeParam)args;
	pCafeTree pcafe = (pCafeTree)param->pcafe;
	double score = 0;
	int skip = 0;
	for ( i = 0 ; i < param->num_params ; i++ )
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
		cafe_set_birthdeath_cache_thread(param);
//		if (param->posterior) {
			score = cafe_get_posterior(param);
//		}
//		else {
//			score = cafe_get_likelihood(param);
//		}
		cafe_free_birthdeath_cache(pcafe);
	}
	char buf[STRING_STEP_SIZE];
	buf[0] = '\0';
	string_pchar_join_double(buf,",", param->num_lambdas, parameters );
	cafe_log(param,"Lambda : %s ", buf, score);
	buf[0] = '\0';
	string_pchar_join_double(buf,",", param->num_mus-param->eqbg, parameters+param->num_lambdas );
	cafe_log(param,"Mu : %s & Score: %f\n", buf, score);
	cafe_log(param, ".");
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
		cafe_set_birthdeath_cache_thread(param);
        //if (param->posterior) {
            score = cafe_get_posterior(param);
        //}
        //else {
        //	score = cafe_get_likelihood(param);
        //}
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
			__cafe_randomize_cluster_parameters( param, param->num_lambdas, param->num_mus, param->k);
            //__cafe_scaleup_cluster_parameters( param, param->num_lambdas, param->num_mus, param->k);
		}
		
		param->pcafe->rootfamilysizes[0] = param->rootfamily_sizes[0];
		param->pcafe->rootfamilysizes[1] = param->rootfamily_sizes[1];
		param->pcafe->familysizes[0] = param->family_sizes[0];
		param->pcafe->familysizes[1] = param->family_sizes[1];
		param->pcafe->rfsize = param->rootfamily_sizes[1] - param->rootfamily_sizes[0] + 1;
		
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
		fminsearch_min(pfm, param->parameters);
		double *re = fminsearch_get_minX(pfm);
		for ( i = 0 ; i < param->num_params ; i++ ) param->parameters[i] = re[i];
        
        
        //__cafe_scaledown_cluster_parameters( param, param->num_lambdas, param->num_mus, param->k);
        

		double current_p = param->parameters[(lambda_len)*(k-param->fixcluster0)];
		double prev_p;
		if (k>0) 
		{
			do {
				double* sumofweights = (double*) memory_new(param->k, sizeof(double));
				for ( i = 0 ; i < param->pfamily->flist->size ; i++ ) {
					for (j = 0; j<k; j++) {
						sumofweights[j] += param->p_z_membership[i][j];
					}
				}
				for (j = 0; j<k-1; j++) {
					param->parameters[(lambda_len)*(k-param->fixcluster0)+j] = sumofweights[j]/param->pfamily->flist->size;
				}
				memory_free(sumofweights);

				fminsearch_min(pfm, param->parameters);	
				
				double *re = fminsearch_get_minX(pfm);
				for ( i = 0 ; i < param->num_params ; i++ ) param->parameters[i] = re[i];
				
				prev_p = current_p;
				current_p = param->parameters[(lambda_len)*(k-param->fixcluster0)];
			} while (current_p - prev_p > pfm->tolx);
		}
		
		cafe_log(param, "\n");
		cafe_log(param,"Lambda Search Result: %d\n", pfm->iters );
		if (k > 0) {
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			if (param->fixcluster0) {
				strncat(buf, "0,", 2);
				string_pchar_join_double(buf,",", param->num_lambdas*(param->k-param->fixcluster0), param->parameters );
			}
			else {
				string_pchar_join_double(buf,",", param->num_lambdas*param->k, param->parameters );
			}
			cafe_log(param,"Lambda : %s\n", buf);
			buf[0] = '\0';
			if (param->k > 0) {
				string_pchar_join_double(buf,",", param->k, param->k_weights );
				cafe_log(param, "p : %s\n", buf);
				cafe_log(param, "p0 : %f\n", param->parameters[param->num_lambdas*(param->k-param->fixcluster0)+0]);
			}
			cafe_log(param, "Score: %f\n", *pfm->fv);
		}
		else {
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			string_pchar_join_double(buf,",", param->num_lambdas, param->parameters );
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
		
		param->pcafe->rootfamilysizes[0] = param->rootfamily_sizes[0];
		param->pcafe->rootfamilysizes[1] = param->rootfamily_sizes[1];
		param->pcafe->familysizes[0] = param->family_sizes[0];
		param->pcafe->familysizes[1] = param->family_sizes[1];
		param->pcafe->rfsize = param->rootfamily_sizes[1] - param->rootfamily_sizes[0] + 1;
		
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
	return param->parameters;
}



double* cafe_best_lambda_mu_by_fminsearch(pCafeParam param, int lambda_len, int mu_len, int k )
{
	int i;
	int max_runs = 10;
	double* scores = memory_new(max_runs, sizeof(double));
	int converged = 0;
	int runs = 0;
	
	do
	{
		if ( param->num_params > 0 )
		{
			__cafe_randomize_cluster_parameters( param, param->num_lambdas, param->num_mus, param->k);
            //__cafe_scaleup_cluster_parameters( param, param->num_lambdas, param->num_mus, param->k);
		}
		
		param->pcafe->rootfamilysizes[0] = param->rootfamily_sizes[0];
		param->pcafe->rootfamilysizes[1] = param->rootfamily_sizes[1];
		param->pcafe->familysizes[0] = param->family_sizes[0];
		param->pcafe->familysizes[1] = param->family_sizes[1];
		param->pcafe->rfsize = param->rootfamily_sizes[1] - param->rootfamily_sizes[0] + 1;
		
		pFMinSearch pfm;
		if (k > 0) {
			pfm = fminsearch_new_with_eq(__cafe_cluster_lambda_mu_search, param->num_params, param);
		}
		else {
			pfm = fminsearch_new_with_eq(__cafe_best_lambda_mu_search, param->num_params, param);
		}
		pfm->tolx = 1e-6;
		pfm->tolf = 1e-6;
		fminsearch_min(pfm, param->parameters);
		double *re = fminsearch_get_minX(pfm);
		for ( i = 0 ; i < param->num_params ; i++ ) param->parameters[i] = re[i];
        
        //__cafe_scaledown_cluster_parameters( param, param->num_lambdas, param->num_mus, param->k);
		
		cafe_log(param, "\n");
		cafe_log(param,"Lambda Search Result: %d\n", pfm->iters );
		// print
		if (k>0) {
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			for( i=0; i<param->num_lambdas; i++) {
				if (param->fixcluster0) {
					strncat(buf, "0,", 2);
					string_pchar_join_double(buf,",", (param->k-param->fixcluster0),  &param->parameters[i*(param->k-param->fixcluster0)] );
				}
				else {
					string_pchar_join_double(buf,",", param->k, &param->parameters[i*param->k] );
				}
				cafe_log(param,"Lambda branch %d: %s\n", i, buf);
				buf[0] = '\0';
			}
			for (i=0; i<param->num_mus-param->eqbg; i++) {
				if (param->fixcluster0) {
					strncat(buf, "0,", 2);
					string_pchar_join_double(buf,",", (param->k-param->fixcluster0),  &param->parameters[param->num_lambdas*(param->k-param->fixcluster0)+i*(param->k-param->fixcluster0)] );
				}
				else {
					string_pchar_join_double(buf,",", param->k, &param->parameters[param->num_lambdas*param->k+i*param->k]);
				}
				cafe_log(param,"Mu branch %d: %s \n", i, buf);
				buf[0] = '\0';
			}
			if (param->k > 0) {
				string_pchar_join_double(buf,",", param->k, param->k_weights );
				cafe_log(param, "p : %s\n", buf);
				cafe_log(param, "p0 : %f\n", param->parameters[param->num_lambdas*(param->k-param->fixcluster0)+(param->num_mus-param->eqbg)*(param->k-param->fixcluster0)+0]);
			}
			cafe_log(param, "Score: %f\n", *pfm->fv);
		}
		else {
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			string_pchar_join_double(buf,",", param->num_lambdas, param->parameters );
			cafe_log(param,"Lambda : %s ", buf, *pfm->fv);
			buf[0] = '\0';
			string_pchar_join_double(buf,",", param->num_mus-param->eqbg, param->parameters+param->num_lambdas );
			cafe_log(param,"Mu : %s & Score: %f\n", buf, *pfm->fv);		
		}
		if (runs > 0) {
			double minscore = __min(scores, runs);
			if (abs(minscore - (*pfm->fv)) < 10*pfm->tolf) {
				converged = 1;
			}
		}
		scores[runs] = *pfm->fv;
		fminsearch_free(pfm);
		
		param->pcafe->rootfamilysizes[0] = param->rootfamily_sizes[0];
		param->pcafe->rootfamilysizes[1] = param->rootfamily_sizes[1];
		param->pcafe->familysizes[0] = param->family_sizes[0];
		param->pcafe->familysizes[1] = param->family_sizes[1];
		param->pcafe->rfsize = param->rootfamily_sizes[1] - param->rootfamily_sizes[0] + 1;
		
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
	return param->parameters;
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
		cafe_set_birthdeath_cache_thread(param);
		double* likelihood = cafe_tree_likelihood(pcafe);
		score = log(__max(likelihood,pcafe->rfsize));
		cafe_free_birthdeath_cache(pcafe);
		pcafe->pbdc_array = NULL;
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

	int rootfamilysizes[2] = { param->rootfamily_sizes[0],  param->rootfamily_sizes[1] };
	int familysizes[2] = { param->family_sizes[0],  param->family_sizes[1] };
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
			pitem->pbdc_array = pref->pbdc_array;
			param->param_set_func(param,pitem->lambda);

			cafe_log(param,"%s: Lambda Search Result of %d/%d in %d iteration \n", pitem->id, i+1, fsize, pfm->iters );
			pString pstr = cafe_tree_string_with_familysize_lambda(param->pcafe);
			cafe_log(param,"%s: %s\n", pitem->id, pstr->buf );
			string_free(pstr);
			continue;
		}

		cafe_family_set_size_with_family_forced(param->pfamily,i,param->pcafe);

		param->rootfamily_sizes[0] = param->pcafe->rootfamilysizes[0];
		param->rootfamily_sizes[1] = param->pcafe->rootfamilysizes[1];
		param->family_sizes[0] = param->pcafe->familysizes[0];
		param->family_sizes[1] = param->pcafe->familysizes[1];

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
//		cafe_set_birthdeath_cache_thread(param);
//		pitem->pbdc_array = param->pcafe->pbdc_array;

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

	param->pcafe->rootfamilysizes[0] = rootfamilysizes[0];
	param->pcafe->rootfamilysizes[1] = rootfamilysizes[1];
	param->pcafe->familysizes[0] = familysizes[0];
	param->pcafe->familysizes[1] = familysizes[1];
	param->pcafe->rfsize = param->rootfamily_sizes[1] - param->rootfamily_sizes[0] + 1;
	param->rootfamily_sizes[0] = rootfamilysizes[0];
	param->rootfamily_sizes[1] = rootfamilysizes[1];
	param->family_sizes[0] = familysizes[0];
	param->family_sizes[1] = familysizes[1];
	memory_free(param->lambda);
	param->lambda = old_lambda;
	return param->lambda;
}

/**************************************************************************
 * Contidional Distribution
**************************************************************************/

typedef struct
{
	pCafeParam cafeparam;
	int range[2];
	pArrayList pCD;
}CDParam;
typedef CDParam* pCDParam;


void* __cafe_conditional_distribution_thread_func(void* ptr)
{
	pCDParam param = (pCDParam)ptr;	
	pCafeParam cafeparam = param->cafeparam;
	pCafeTree pcafe = cafe_tree_copy(cafeparam->pcafe);
	pcafe->pbdc_array = cafeparam->pcafe->pbdc_array;
#ifdef __DEBUG__
	printf("CD: %d ~ %d\n", param->range[0], param->range[1]);
#endif
	param->pCD = cafe_tree_conditional_distribution(pcafe, param->range, cafeparam->num_random_samples);
	cafe_tree_free(pcafe);
	return (NULL);
}

pArrayList cafe_conditional_distribution(pCafeParam param)
{
	int numthreads = param->num_threads;
	int threadstep = param->pcafe->rfsize/numthreads;
	if ( threadstep == 0 )
	{
		numthreads = param->pcafe->rfsize;
	}
	else
	{
		threadstep--;
	}

	pCDParam ptparam = (pCDParam)memory_new(numthreads,sizeof(CDParam));
	int i, r = param->rootfamily_sizes[0] ;
	for ( i = 0 ; i < numthreads; i++, r+=threadstep+1 )
	{
		ptparam[i].cafeparam = param;
		ptparam[i].range[0]= r;
		ptparam[i].range[1]= r + threadstep;
	}
	ptparam[numthreads-1].range[1] = param->rootfamily_sizes[1];
	thread_run(numthreads, __cafe_conditional_distribution_thread_func, ptparam, sizeof(CDParam));
	pArrayList cdlist = ptparam[0].pCD;
	for( i = 1 ; i < numthreads ; i++ )
	{
		for ( r = 0 ; r < ptparam[i].pCD->size ; r++ )
		{
			arraylist_add(cdlist, ptparam[i].pCD->array[r]);
		}
		arraylist_free(ptparam[i].pCD, NULL);
	}
	memory_free(ptparam);
	ptparam = NULL;
	return cdlist;
}

/**************************************************************************
 * Viterbi
**************************************************************************/

typedef struct
{
	pCafeParam cafeparam;
	pArrayList pCD;
	int from;
}ViterbiParam;

typedef ViterbiParam*  pViterbiParam;

pthread_mutex_t mutex_cafe_viterbi = PTHREAD_MUTEX_INITIALIZER;

void* __cafe_viterbi_thread_func(void* ptr)
{
	int i, j, k, m;
	pViterbiParam pv = (pViterbiParam)ptr;
	pCafeParam param = (pCafeParam)pv->cafeparam;
	pArrayList pCD = pv->pCD;
	pCafeTree pcafe = cafe_tree_copy(param->pcafe);
	pcafe->pbdc_array = param->pcafe->pbdc_array;
	pTree ptree = (pTree)pcafe;
	int nnodes = (ptree->nlist->size-1)/2;
	int fsize = param->pfamily->flist->size;
	double* cP = (double*)memory_new( pcafe->rfsize, sizeof(double));
#ifdef __DEBUG__
	printf("VITERBI: from %d\n", pv->from );
#endif
	for ( i = pv->from; i < fsize ; i+=param->num_threads )
	{
		cafe_family_set_size_with_family_forced(param->pfamily,i, pcafe);
//		cafe_family_set_size(param->pfamily,i, pcafe);
		
		cafe_tree_p_values(pcafe, cP, pCD, param->num_random_samples);
		param->viterbi.maximumPvalues[i] = __max(cP,pcafe->rfsize);
		cafe_tree_viterbi(pcafe);
		/* check family size for all nodes first */
		for ( j = 0 ; j < nnodes ; j++ )
		{
			pCafeNode pcnode = (pCafeNode)ptree->nlist->array[j];
			if(pcnode->familysize>10000) {
				fprintf(stderr,"ERROR: FamilySize larger than bd array size Something wrong\n"); 
				exit(-1);
			}
		}
		/* end check family size for all nodes first */

		for ( j = 0 ; j < nnodes ; j++ )
		{
			pCafeNode pcnode = (pCafeNode)ptree->nlist->array[2*j+1];
			param->viterbi.viterbiNodeFamilysizes[j][i] = pcnode->familysize;
			pCafeNode child[2] = { (pCafeNode)((pTreeNode)pcnode)->children->head->data, 
							       (pCafeNode)((pTreeNode)pcnode)->children->tail->data };
			for ( k = 0 ; k < 2 ; k++ )
			{
				m = j*2 + k;
				if ( child[k]->familysize > pcnode->familysize ) param->viterbi.expandRemainDecrease[0][m]++;
				else if ( child[k]->familysize == pcnode->familysize ) param->viterbi.expandRemainDecrease[1][m]++;
				else param->viterbi.expandRemainDecrease[2][m]++;
pthread_mutex_lock( &mutex_cafe_viterbi );
				param->viterbi.averageExpansion[m] += child[k]->familysize - pcnode->familysize;
pthread_mutex_unlock( &mutex_cafe_viterbi );
			}
		}

		if ( param->viterbi.maximumPvalues[i] > param->pvalue )
		{
			for ( j = 0 ; j < ptree->nlist->size-1 ; j++ )
			{
				param->viterbi.viterbiPvalues[j][i] = -1;
			}
			continue;
		}

		for ( j = 0 ; j < nnodes ; j++ )
		{
			pCafeNode pcnode = (pCafeNode)ptree->nlist->array[2*j+1];
			pCafeNode child[2] = { (pCafeNode)((pTreeNode)pcnode)->children->head->data, 
							       (pCafeNode)((pTreeNode)pcnode)->children->tail->data };
			for ( k = 0 ; k < 2 ; k++ )
			{
				double p = child[k]->bd[pcnode->familysize][child[k]->familysize];
				double** pbdc = child[k]->bd;
				int n = 2 * j + k;
				for ( m = 0 ; m <= pcafe->familysizes[1]; m++ )
				{
					if ( pbdc[pcnode->familysize][m] == p )
					{
						param->viterbi.viterbiPvalues[n][i] += pbdc[pcnode->familysize][m]/2.0;
					}
					else if ( pbdc[pcnode->familysize][m] < p )
					{
						param->viterbi.viterbiPvalues[n][i] += pbdc[pcnode->familysize][m];
					}
				}
			}
		}
	}
	memory_free(cP);
	cP = NULL;
	cafe_tree_free(pcafe);
	return (NULL);
}

pArrayList cafe_viterbi(pCafeParam param, pArrayList pCD)
{
	cafe_log(param,"Running Viterbi algorithm....\n");

	pViterbiParam ptparam = (pViterbiParam)memory_new(param->num_threads,sizeof(ViterbiParam));
	pTree ptree = (pTree)param->pcafe;

	if ( pCD == NULL )
	{
		param->param_set_func(param,param->parameters);
		cafe_set_birthdeath_cache_thread(param);
		pCD = cafe_conditional_distribution(param);
		//cafe_free_birthdeath_cache(param->pcafe);
	}

	int nrows = param->pfamily->flist->size;
	int nnodes = ptree->nlist->size - 1;

	param->viterbi.viterbiPvalues = (double**)memory_new_2dim(nnodes,nrows,sizeof(double));
	param->viterbi.expandRemainDecrease = (int**)memory_new_2dim(3,nnodes,sizeof(int));
	param->viterbi.viterbiNodeFamilysizes = (int**)memory_new_2dim(nnodes,nrows,sizeof(int));
	param->viterbi.maximumPvalues = (double*)memory_new(nrows, sizeof(double));
	param->viterbi.averageExpansion = (double*)memory_new(nnodes, sizeof(double));

	int i;
	for ( i = 0 ; i < param->num_threads; i++ )
	{
		ptparam[i].cafeparam = param;
		ptparam[i].from = i;
		ptparam[i].pCD = pCD;
	}
	thread_run(param->num_threads, __cafe_viterbi_thread_func, ptparam, sizeof(ViterbiParam));

	for ( i = 0 ; i < ptree->nlist->size - 1; i++ )
	{
		param->viterbi.averageExpansion[i] /= param->pfamily->flist->size;
	}

	memory_free(ptparam);	
	ptparam = NULL;				
	return pCD;
}

void cafe_viterbi_print(pCafeParam param)
{
	int i, j;
	int size = param->pfamily->flist->size;
	pCafeTree pcafe = param->pcafe;
	pTree ptree = (pTree)pcafe;
	for ( i = 0 ; i < size ; i++ )
	{
		cafe_family_set_size(param->pfamily,i, pcafe);
		for ( j = 1 ; j < ptree->nlist->size ; j+=2 )
		{
			pCafeNode pcnode = (pCafeNode)ptree->nlist->array[j];
			pcnode->familysize = param->viterbi.viterbiNodeFamilysizes[j/2][i];
		}
		cafe_tree_string_print(pcafe);
	}
}

/**************************************************************************
 * BranchCutting
**************************************************************************/
typedef struct
{
	pCafeParam cafeparam;
	pArrayList** pCDSs;				
	int range[2];
}BranchCuttingParam;

typedef BranchCuttingParam* pBranchCuttingParam;

void* __cafe_branch_cutting_thread_func(void* ptr)
{
	int i, b;
	pBranchCuttingParam pbc = (pBranchCuttingParam)ptr;
	pCafeParam param = pbc->cafeparam;	
	pArrayList** pCDSs = pbc->pCDSs;

#ifdef __DEBUG__
	printf("Branch cutting : %d ~ %d\n", pbc->range[0], pbc->range[1] -1 );
#endif

	pTree ptree = (pTree)param->pcafe;
	int nnodes = ptree->nlist->size;
	double* p1 = (double*)memory_new( param->pcafe->rfsize, sizeof(double));
	double** p2 = (double**)memory_new_2dim( param->pcafe->rfsize, param->pcafe->rfsize, sizeof(double));

	for ( b = 0 ; b < nnodes ; b++ )
	{
		if ( tree_is_root( ptree,(pTreeNode)ptree->nlist->array[b] ) ) 
		{
			continue;
		}
		
		pCafeTree pcafe = cafe_tree_copy(param->pcafe);
		pcafe->pbdc_array = param->pcafe->pbdc_array;
		pCafeTree psub =  cafe_tree_split(pcafe,b);
		psub->pbdc_array = param->pcafe->pbdc_array;

		for ( i = pbc->range[0] ; i < pbc->range[1] ; i++ )
		{
			pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
			if ( pitem->ref >= 0 && pitem->ref != i ) continue;
			if ( param->viterbi.maximumPvalues[i] > param->pvalue )
			{
				param->cutPvalues[b][i] = -1;
				continue;
			}
			if ( tree_is_leaf( psub->super.root ) || tree_is_leaf(pcafe->super.root) )
			{
				pCafeTree pct = tree_is_leaf ( psub->super.root ) ? pcafe : psub;
				cafe_family_set_size_for_split(param->pfamily,i, pct);
				cafe_tree_p_values( pct,p1,pCDSs[b][0], param->num_random_samples  );
				param->cutPvalues[b][i] = __max(p1, pcafe->rfsize );
			}
			else
			{
				cafe_family_set_size_for_split(param->pfamily,i, pcafe);
				cafe_family_set_size_for_split(param->pfamily,i, psub);
				cafe_tree_p_values_of_two_trees(pcafe,psub, p2,
										pCDSs[b][0], pCDSs[b][1], param->num_random_samples/10 );
				param->cutPvalues[b][i] = 0;
				int m,n;
				for ( m = 0 ; m < pcafe->rfsize ; m++ )
				{
					for( n = 0 ; n < pcafe->rfsize ; n++ )
					{
						if ( p2[m][n] > param->cutPvalues[b][i] )
						{
							param->cutPvalues[b][i] = p2[m][n];
						}
					}
				}
			}
		}
		cafe_tree_free(pcafe);	
		cafe_tree_free(psub);
	}
	memory_free(p1);
	p1 = NULL;
	memory_free_2dim((void**)p2,param->pcafe->rfsize,param->pcafe->rfsize,NULL);
	return (NULL);
}

void cafe_branch_cutting(pCafeParam param)
{
	cafe_log(param,"Running Branch Cutting....\n");

	pTree ptree = (pTree)param->pcafe;
	int i,b,j;
	int nnodes = ptree->nlist->size;
	pArrayList** pCDSs = (pArrayList**)memory_new_2dim(nnodes,2,sizeof(pArrayList));

	for ( b = 0 ; b < nnodes  ; b++ )
	{
		if ( tree_is_root( ptree,(pTreeNode)ptree->nlist->array[b] ) ) 
		{
			pCDSs[b][0] = NULL;
			pCDSs[b][1] = NULL;
			continue;
		}
		pCafeTree pcafe = cafe_tree_copy(param->pcafe);
		pcafe->pbdc_array = param->pcafe->pbdc_array;
		pCafeTree psub =  cafe_tree_split(pcafe,b);
		psub->pbdc_array = param->pcafe->pbdc_array;
		
		cafe_log(param,">> %d  --------------------\n", b ); 
		pString pstr = cafe_tree_string(pcafe);
		cafe_log(param,"%s\n", pstr->buf );
		string_free(pstr);
		pstr = cafe_tree_string(psub);
		cafe_log(param,"%s\n", pstr->buf );
		string_free(pstr);

		pCafeTree porig = param->pcafe;
		if ( tree_is_leaf( psub->super.root ) )
		{
			param->pcafe = pcafe;
			pCDSs[b][0] = cafe_conditional_distribution(param);
			pCDSs[b][1] = NULL;
		}
		else if ( tree_is_leaf ( pcafe->super.root) )
		{
			param->pcafe = psub;
			pCDSs[b][0] = cafe_conditional_distribution(param);
			pCDSs[b][1] = NULL;
		}
		else
		{
			int orig_r = param->num_random_samples;
			param->num_random_samples /= 10;
			param->pcafe = pcafe;
			pCDSs[b][0] = cafe_conditional_distribution(param);
			param->pcafe = psub;
			pCDSs[b][1] = cafe_conditional_distribution(param);
			param->num_random_samples = orig_r;
		}
		param->pcafe = porig;

		cafe_tree_free(pcafe);	
		cafe_tree_free(psub);
	}

	int threadstep = param->pfamily->flist->size/param->num_threads;
	pBranchCuttingParam ptparam = (pBranchCuttingParam)memory_new(param->num_threads,sizeof(BranchCuttingParam));

	int nrows = param->pfamily->flist->size;
	param->cutPvalues = (double**)memory_new_2dim(nnodes,nrows,sizeof(double));
	int rid = ptree->root->id;
	for ( i = 0 ; i < nrows; i++ )
	{
		param->cutPvalues[rid][i] = -1;
	}

	int r = 0;
	for ( i = 0 ; i < param->num_threads; i++, r+=threadstep )
	{
		ptparam[i].cafeparam = param;
		ptparam[i].pCDSs = pCDSs;
		ptparam[i].range[0] = r;
		ptparam[i].range[1] = r + threadstep;
	}
	ptparam[param->num_threads-1].range[1]= param->pfamily->flist->size;
	thread_run(param->num_threads, __cafe_branch_cutting_thread_func, ptparam, sizeof(BranchCuttingParam));

	for ( i = 0 ; i < nrows ; i++ )
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
		if ( pitem->ref < 0 || pitem->ref == i ) continue;
		for( b = 0 ; b < nnodes ; b++ )
		{
			param->cutPvalues[b][i] = param->cutPvalues[b][pitem->ref];
		}
	}

	memory_free(ptparam);		
	ptparam = NULL;			
	for ( i = 0 ; i < nnodes; i++ )
	{
		for ( j = 0 ; j < 2 ; j++ )
		{
			if ( pCDSs[i][j] ) arraylist_free(pCDSs[i][j], free);
		}
	}
	memory_free(pCDSs);
	pCDSs = NULL;
	cafe_log( param , "Done : Branch Cutting\n" );
}

/**************************************************************************
 * Likelihood ratio test
**************************************************************************/
typedef struct
{
	pCafeParam cafeparam;
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
	pcafe->pbdc_array = param->pcafe->pbdc_array;
	int nnodes = ptree->nlist->size;
	int old_bl;
	int fsize = param->pfamily->flist->size;
#ifdef __DEBUG__
	printf("Likelihood ratio test: %d\n", plrt->from );
#endif
	for ( i = plrt->from ; i < fsize ; i += param->num_threads )
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
		if ( pitem->ref >= 0 &&  pitem->ref != i ) continue;
		if ( param->viterbi.maximumPvalues[i] > param->pvalue )
		{
			for( b = 0 ; b < nnodes ; b++ ) param->likelihoodRatios[b][i] = -1;
			continue;
		}
		cafe_family_set_size(param->pfamily,i, pcafe);
		double maxlh = __max(cafe_tree_likelihood(pcafe), param->pcafe->rfsize );
		for( b = 0 ; b < nnodes ; b++ )
		{
			pPhylogenyNode pnode = (pPhylogenyNode)pcafe->super.nlist->array[b];
			if ( tree_is_root( (pTree)pcafe ,(pTreeNode)pnode) ) 
			{
				param->likelihoodRatios[b][i] = -1;
				continue;
			}
			old_bl = pnode->branchlength;
			double** old_bd = ((pCafeNode)pnode)->bd;
			double prevlh = -1;
			double nextlh = maxlh;
			while( prevlh < nextlh )
			{
				prevlh = nextlh;
				pnode->branchlength += rint(pnode->branchlength * 0.15);
pthread_mutex_lock( &mutex_cafe_likelihood );
				((pCafeNode)pnode)->bd = birthdeath_cache_get_matrix(pcafe->pbdc_array, pnode->branchlength, ((pCafeNode)pnode)->lambda,  ((pCafeNode)pnode)->mu );
pthread_mutex_unlock( &mutex_cafe_likelihood);
				nextlh = __max(cafe_tree_likelihood(pcafe), param->pcafe->rfsize );
			}
			param->likelihoodRatios[b][i] = (prevlh == maxlh) ? 1 : 1 - chi2cdf( 2*(log(prevlh) - log(maxlh)), 1);
		//	param->likelihoodRatios[b][i] = (prevlh == maxlh) ? 1 : prevlh / maxlh ;
			pnode->branchlength = old_bl;
			((pCafeNode)pnode)->bd = old_bd;
		}
	}
	cafe_tree_free(pcafe);
	return (NULL);
}

void cafe_likelihood_ratio_test(pCafeParam param)
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

/**************************************************************************
 * Likelihood ratio test with more than one lambda
**************************************************************************/

typedef struct
{
	pCafeParam cafeparam;
	int from;
	int* lambda;
	double* pvalues;
	param_func lfunc;
	pTree lambda_tree;
	int    num_lambdas;
}LRT2LParam;
typedef LRT2LParam* pLRT2LParam;

double** lambda_cache;
pBirthDeathCacheArray* PBDC;

pthread_mutex_t mutex_cafe_lh2 = PTHREAD_MUTEX_INITIALIZER;

double __cafe_lhr_get_likelihood_for_diff_lambdas(pCafeParam param, int idx, int t)
{
	param->branchlength_update_func(param,&t);
pthread_mutex_lock( &mutex_cafe_lh2 );
	if ( lambda_cache[t] == NULL )
	{
		param->lambda = NULL;
		cafe_best_lambda_by_fminsearch(param,param->num_lambdas, 0);
		lambda_cache[t] = param->lambda;
		cafe_set_birthdeath_cache(param);
		PBDC[t] = param->pcafe->pbdc_array;
	}
	else
	{
		memcpy( param->lambda, lambda_cache[t], sizeof(double)*param->num_lambdas );
		param->param_set_func(param,param->lambda );
		param->pcafe->pbdc_array = PBDC[t];
		cafe_tree_set_birthdeath(param->pcafe);
	}
pthread_mutex_unlock( &mutex_cafe_lh2 );
	int i;
	cafe_family_set_size(param->pfamily,idx, param->pcafe);
	double mlh = __max(cafe_tree_likelihood(param->pcafe), param->pcafe->rfsize );
	pTree ptree = (pTree)param->pcafe;
	pArrayList nlist = ptree->nlist;
	for ( i = 0 ; i < nlist->size ; i++ )
	{
		pPhylogenyNode pnode = (pPhylogenyNode) nlist->array[i];
		pnode->branchlength = param->old_branchlength[i];
	}
	return mlh;
}

void* __cafe_lhr_for_diff_lambdas_thread(void* ptr)
{
	int i, j;
	pLRT2LParam plrt = (pLRT2LParam)ptr;
	pCafeParam param = plrt->cafeparam;
	pCafeTree pcafe = cafe_tree_copy(param->pcafe);
	pcafe->pbdc_array = param->pcafe->pbdc_array;

	pCafeParam cpy_param = cafe_copy_parameters(param);	
	cpy_param->num_lambdas = plrt->num_lambdas;
	cpy_param->lambda_tree = plrt->lambda_tree;
	cpy_param->param_set_func = plrt->lfunc;

	int fsize = param->pfamily->flist->size;
	for( i = plrt->from; i < fsize ; i+= param->num_threads )
	{ 
		pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
		if ( pitem->ref >= 0 && pitem->ref != i ) continue;
		cafe_family_set_size(param->pfamily,i, pcafe);
		double maxlh1 = __max(cafe_tree_likelihood(pcafe), pcafe->rfsize );
		double prev = -1;
		double next = __cafe_lhr_get_likelihood_for_diff_lambdas(cpy_param,i,0);
		for ( j = 1 ; prev < next ; j++ )
		{
			prev = next;
			next = __cafe_lhr_get_likelihood_for_diff_lambdas(cpy_param,i,j);
		}
		plrt->pvalues[i] = (prev == maxlh1 ) ? 1 :  2*(log(prev) - log(maxlh1));
		plrt->lambda[i] = j - 2;
	}
	cafe_free_copy_parameters(cpy_param);
	cafe_tree_free(pcafe);
	return (NULL);	
}

void cafe_lhr_report(pCafeParam param, double* pvalues, int* plambda )
{
	int i;
	int fsize = param->pfamily->flist->size;
	for ( i = 0 ; i < fsize ; i++ )
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
		cafe_family_set_size(param->pfamily,i, param->pcafe);
		pString pstr = cafe_tree_string(param->pcafe);
		fprintf(param->fout, "%s\t%s\t", pitem->id, pstr->buf );		
		double* l = lambda_cache[plambda[i]];
		fprintf(param->fout, "(%d, %lf,%lf)\t%g\t%f\n", plambda[i], l[0],l[1], pvalues[i], pvalues[1] == 1 ? 1 : 1-chi2cdf(pvalues[i],1));
		string_free(pstr);
	}
}

void cafe_lhr_for_diff_lambdas(pCafeParam param, pTree lambda_tree2, int num_lambdas, param_func lfunc )
{
	cafe_log(param,"Running Likelihood Ratio Test 2....\n");
	int i;
	lambda_cache = (double**)memory_new(100,sizeof(double*));
	PBDC = (pBirthDeathCacheArray*)memory_new(100,sizeof(pBirthDeathCacheArray));
	for ( i = 0 ; i < 100 ; i++ )
	{
		lambda_cache[i] = NULL;
		PBDC[i] = NULL;
	}
	int old_numthreads = param->num_threads;
	param->num_threads = 1;
	int nrows = param->pfamily->flist->size;
	double* pvalues = (double*)memory_new(nrows,sizeof(double));
	int* plambda = (int*)memory_new(nrows, sizeof(int));

	pLRT2LParam ptparam = (pLRT2LParam) memory_new( param->num_threads, sizeof(LRT2LParam) );
	for( i = 0 ; i < param->num_threads ; i++ )
	{
		ptparam[i].cafeparam = param;
		ptparam[i].from = i;
		ptparam[i].lambda = plambda;
		ptparam[i].pvalues = pvalues;
		ptparam[i].lfunc = lfunc;
		ptparam[i].lambda_tree = lambda_tree2;
		ptparam[i].num_lambdas = num_lambdas;
	}
	param_func old_func = param->param_set_func;
	param->param_set_func = cafe_lambda_set_default;
	thread_run(param->num_threads, __cafe_lhr_for_diff_lambdas_thread, ptparam, sizeof(LRT2LParam));
	param->param_set_func = old_func;

	int fsize = param->pfamily->flist->size;
	for ( i = 0 ; i < fsize ; i++ )
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
		if ( pitem->ref < 0 || pitem->ref == i ) continue;
		pvalues[i] = pvalues[ pitem->ref ];
		plambda[i] = plambda[pitem->ref];
	}

	cafe_lhr_report(param, pvalues, plambda );

	for ( i = 0 ; i < 100 ; i++ )
	{
		if ( lambda_cache[i] ) 
		{
			memory_free(lambda_cache[i]);
			lambda_cache[i] = NULL;
			birthdeath_cache_array_free(PBDC[i]);
		}
	}
	memory_free(pvalues);
	pvalues = NULL;
	memory_free(plambda);
	plambda = NULL;
	param->num_threads = old_numthreads;
	memory_free(lambda_cache);
	lambda_cache = NULL;
	memory_free(PBDC);
	PBDC = NULL;
}

/*************************************************************************
 * Deprecated
 *************************************************************************/

/**************************************************************************
 * Main
**************************************************************************/

void* cafe_run(void* ptr)
{
	pCafeParam param = (pCafeParam)ptr;
	if ( param->lambda == NULL )
	{
		cafe_best_lambda_by_fminsearch(param, param->num_lambdas, 0);
	}
	else
	{
		param->param_set_func(param,param->lambda);
		pString pstr = cafe_tree_string_with_lambda(param->pcafe);
		cafe_log(param, "Lambda Value: %s\n", pstr->buf );
		string_free(pstr);
	}
	cafe_set_birthdeath_cache(param);
	pArrayList pCD = cafe_viterbi(param, NULL);
//	cafe_branch_cutting(param);
	cafe_likelihood_ratio_test(param);
	cafe_free_birthdeath_cache(param->pcafe);
	cafe_report(param, CAFE_REPORT_TEXT);
	arraylist_free(pCD, free);
	return (NULL);
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
	//param->branchlengths_sorted = (int*)memory_new(psrc->num_branches,sizeof(int));
	//memcpy( param->branchlengths_sorted, psrc->branchlengths_sorted, psrc->num_branches * sizeof(int));

	param->viterbi.viterbiPvalues = NULL;
	param->viterbi.expandRemainDecrease = NULL;
	param->viterbi.viterbiNodeFamilysizes = NULL;
	param->viterbi.maximumPvalues = NULL;
	param->viterbi.averageExpansion = NULL;
	param->cutPvalues = NULL;

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
