/** @file cafe_shell.c
* @brief  Functions corresponding to the commands available in Cafe.
*
* Command list is found in the #cafe_cmd function
*/

#include<mathfunc.h>
#include<ctype.h>
#include <assert.h>
#include <float.h>
#include<stdarg.h>
#include<stdio.h>
#include "cafe.h"
#include "cafe_shell.h"

/**
* \brief Holds the global program state that user commands act on.
*
*/
pCafeParam cafe_param;

pTree tmp_lambda_tree;

#ifndef STDERR_IF
	#define STDERR_IF(a,b)	if ( a ) { fprintf(stderr,b); return -1; }
#endif

void cafe_shell_set_lambda(pCafeParam param, double* lambda);
void cafe_shell_set_lambda_mu(pCafeParam param, double* parameters);

void cafe_shell_set_lambdas(pCafeParam param, double* lambda)
{
    if (param->optimizer_init_type == LAMBDA_ONLY)
        cafe_shell_set_lambda(param, lambda);

    if (param->optimizer_init_type == LAMBDA_MU)
        cafe_shell_set_lambda_mu(param, lambda);
}

void reset_k_likelihoods(pCafeNode pcnode, int k, int num_factors)
{
	if (pcnode->k_likelihoods) { memory_free(pcnode->k_likelihoods); pcnode->k_likelihoods = NULL; }
	pcnode->k_likelihoods = (double**)memory_new_2dim(k, num_factors, sizeof(double));
}

void set_birth_death_probabilities(struct probabilities *probs, int num_lambdas, int first_mu, int fix_cluster, double* parameters)
{
	if (num_lambdas < 1)
	{
		probs->lambda = parameters[0];
		probs->mu = parameters[first_mu];
	}
	else
	{
		probs->lambda = -1;
		probs->mu = -1;
		free_probabilities(probs);
		probs->param_lambdas = (double*)memory_new(num_lambdas, sizeof(double));
		if (!fix_cluster) {
			memcpy(&probs->param_lambdas[0], &parameters[0], (num_lambdas)*sizeof(double));
		}
		else {
			probs->param_lambdas[0] = 0;
			memcpy(&probs->param_lambdas[1], &parameters[0], (num_lambdas - fix_cluster)*sizeof(double));
		}

		probs->param_mus = (double*)memory_new(num_lambdas, sizeof(double));
		if (!fix_cluster) {
			memcpy(&probs->param_mus[0], &parameters[first_mu*num_lambdas], (num_lambdas)*sizeof(double));
		}
		else {
			probs->param_mus[0] = 0;
			memcpy(&probs->param_mus[1], &parameters[first_mu*(num_lambdas - fix_cluster)], (num_lambdas - fix_cluster)*sizeof(double));
		}
	}
}

void set_birth_death_probabilities2(struct probabilities *probs, int num_lambdas, int first_mu, int fix_cluster, int taxa_id, int eqbg, double* parameters)
{
	if (taxa_id < 0)
	{
#ifdef VERBOSE
		fprintf(stderr, "WARNING: No taxa id provided, assuming 0\n");
#endif
		taxa_id = 0;
	}

	if (num_lambdas > 0) {
		probs->lambda = -1;
		probs->mu = -1;

		free_probabilities(probs);
		// set lambdas
		probs->param_lambdas = (double*)memory_new(num_lambdas, sizeof(double));
		if (!fix_cluster) {
			memcpy(&probs->param_lambdas[0], &parameters[taxa_id*num_lambdas], (num_lambdas)*sizeof(double));
		}
		else {
			probs->param_lambdas[0] = 0;
			memcpy(&probs->param_lambdas[1], &parameters[taxa_id*(num_lambdas - 1)], (num_lambdas - 1)*sizeof(double));
		}

		// set mus
		probs->param_mus = (double*)memory_new(num_lambdas, sizeof(double));
		if (eqbg) {
			if (taxa_id == 0) {
				memcpy(probs->param_mus, probs->param_lambdas, (num_lambdas - fix_cluster)*sizeof(double));
			}
			else {
				if (!fix_cluster) {
					memcpy(&probs->param_mus[0], &parameters[(first_mu)*num_lambdas + (taxa_id - eqbg)*num_lambdas], (num_lambdas)*sizeof(double));
				}
				else {
					probs->param_mus[0] = 0;
					memcpy(&probs->param_mus[1], &parameters[(first_mu)*(num_lambdas - 1) + (taxa_id - eqbg)*(num_lambdas - 1)], (num_lambdas - 1)*sizeof(double));
				}
			}
		}
		else {
			if (!fix_cluster) {
				memcpy(&probs->param_mus[0], &parameters[(first_mu)*num_lambdas + taxa_id*num_lambdas], (num_lambdas)*sizeof(double));
			}
			else {
				probs->param_mus[0] = 0;
				memcpy(&probs->param_mus[1], &parameters[(first_mu)*(num_lambdas - 1) + taxa_id*(num_lambdas - 1)], (num_lambdas - 1)*sizeof(double));
			}
		}
	}
	else
	{
		if (eqbg) {
			probs->lambda = parameters[taxa_id];
			if (taxa_id == 0) {
				probs->mu = probs->lambda;
			}
			else {
				probs->mu = parameters[(first_mu)+(taxa_id - eqbg)];
			}
		}
		else {
			probs->lambda = parameters[taxa_id];
			probs->mu = parameters[(first_mu)+taxa_id];
		}

	}
}

void set_birth_death_probabilities4(struct probabilities *probs, int num_lambdas, int fix_cluster, int taxa_id, double* parameters)
{
	if (taxa_id < 0)
	{
#ifdef VERBOSE
		fprintf(stderr, "WARNING: No taxa id provided, assuming 0\n");
#endif
		taxa_id = 0;
	}
	if (num_lambdas > 0) {
		probs->lambda = -1;
		probs->mu = -1;
		free_probabilities(probs);

		probs->param_lambdas = (double*)memory_new(num_lambdas, sizeof(double));
		if (!fix_cluster) {
			memcpy(&probs->param_lambdas[0], &parameters[taxa_id*num_lambdas], (num_lambdas)*sizeof(double));
		}
		else {
			// TODO: Investigate this. Subtracting 1 from num_lambdas does not seem like correct behavior.
			probs->param_lambdas[0] = 0;
			memcpy(&probs->param_lambdas[1], &parameters[taxa_id*(num_lambdas - 1)], (num_lambdas - 1)*sizeof(double));
		}

	}
	else {
		probs->lambda = parameters[taxa_id];
		probs->mu = -1;
	}
}

void initialize_z_membership(pCafeParam param)
{
	if (param->p_z_membership == NULL) {
		param->p_z_membership = (double**)memory_new_2dim(param->pfamily->flist->size, param->num_lambdas*param->parameterized_k_value, sizeof(double));
		// assign based on param->k_weights (prior)
		for (int i = 0; i < param->pfamily->flist->size; i++)
		{
			for (int k = 0; k<param->parameterized_k_value; k++) {
				param->p_z_membership[i][k] = param->k_weights[k];
			}
		}
	}

}

void initialize_k_bd(pCafeParam param, double *parameters)
{
	pArrayList nlist = param->pcafe->super.nlist;
	for (int i = 0; i < nlist->size; i++)
	{
		int taxa_id = 0;
		if (param->lambda_tree != NULL)
			taxa_id = ((pPhylogenyNode)param->lambda_tree->nlist->array[i])->taxaid;

		pCafeNode pcnode = (pCafeNode)nlist->array[i];

		if (param->parameterized_k_value > 0) {
			reset_k_likelihoods(pcnode, param->parameterized_k_value, param->pcafe->size_of_factor);

			if (pcnode->k_bd) { arraylist_free(pcnode->k_bd, NULL); }
			pcnode->k_bd = arraylist_new(param->parameterized_k_value);
		}
		set_birth_death_probabilities4(&pcnode->birth_death_probabilities, param->parameterized_k_value, param->fixcluster0, taxa_id, parameters);
	}
}

void initialize_k_bd2(pCafeParam param, double *parameters)
{
	pArrayList nlist = param->pcafe->super.nlist;
	for (int i = 0; i < nlist->size; i++)
	{
		int reset = 0;
		pCafeNode pcnode = (pCafeNode)nlist->array[i];
		if (param->lambda_tree != NULL)
		{
			int taxa_id = ((pPhylogenyNode)param->lambda_tree->nlist->array[i])->taxaid;
			set_birth_death_probabilities2(&pcnode->birth_death_probabilities, param->parameterized_k_value, param->num_lambdas, param->fixcluster0, taxa_id, param->eqbg, parameters);
			reset = 1;
		}
		else
		{
			set_birth_death_probabilities(&pcnode->birth_death_probabilities, param->parameterized_k_value, param->num_lambdas, param->fixcluster0, parameters);
			if (param->parameterized_k_value > 0)
				reset = 1;
		}

		if (reset != 0)
		{
			reset_k_likelihoods(pcnode, param->parameterized_k_value, param->pcafe->size_of_factor);

			if (pcnode->k_bd) { arraylist_free(pcnode->k_bd, NULL); }
			pcnode->k_bd = arraylist_new(param->parameterized_k_value);
		}
	}

}
void cafe_shell_set_lambda(pCafeParam param, double* parameters)
{
	if (param->input.parameters[0] != parameters[0])
		memcpy(param->input.parameters, parameters, param->num_params*sizeof(double));
	
	// set lambda
	param->lambda = param->input.parameters;
	
	// set k_weights
	if (param->parameterized_k_value > 0) {
		int start = param->num_lambdas*(param->parameterized_k_value - param->fixcluster0);
		input_values_copy_weights(param->k_weights, &param->input, start, param->parameterized_k_value);
		initialize_z_membership(param);
	}
	
	param->pcafe->k = param->parameterized_k_value;
	initialize_k_bd(param, parameters);
}

void cafe_shell_set_lambda_mu(pCafeParam param, double* parameters)
{
	if (param->input.parameters[0] != parameters[0]) {
		memcpy(param->input.parameters, parameters, param->num_params*sizeof(double));
	}
	// set lambda and mu
	param->lambda = param->input.parameters;
	if (param->parameterized_k_value > 0) {
		param->mu = &(param->input.parameters[param->num_lambdas*(param->parameterized_k_value-param->fixcluster0)]);
	}
	else {
		param->mu = &(param->input.parameters[param->num_lambdas]);
	}

	// set k_weights
	if (param->parameterized_k_value > 0) {
		int start = param->num_lambdas*(param->parameterized_k_value - param->fixcluster0) + (param->num_mus - param->eqbg)*(param->parameterized_k_value - param->fixcluster0);
		input_values_copy_weights(param->k_weights, &param->input, start, param->parameterized_k_value);
		initialize_z_membership(param);
	}
	
	param->pcafe->k = param->parameterized_k_value;
	initialize_k_bd2(param, parameters);
}

void cafe_shell_set_branchlength(pCafeParam param, int max_family_size)
{
	int i;
	char buf[STRING_STEP_SIZE];

	pArrayList nlist = param->pcafe->super.nlist;
	for ( i = 0; i < nlist->size ; i++ )
	{
		pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
		if ( tree_is_root( (pTree)param->pcafe, (pTreeNode)pnode) ) continue;
		printf("%d[%d]: ", i, (int)pnode->branchlength );
		if (fgets(buf,STRING_STEP_SIZE,stdin) == NULL)
			fprintf(stderr, "Failed to read input\n");

		size_t len = strlen(buf);
		buf[--len] = '\0';
		if ( len != 0 )
		{
			int size = -1;
			sscanf( buf, "%d", &size );
			if ( size > 0 )
			{ 
				pnode->branchlength = size;
			}
			else
			{
				fprintf(stderr,"ERROR: the branch length of node %d is not changed\n", i);
			}
		}
	}
	
    cafe_tree_set_birthdeath(param->pcafe, max_family_size);
}


void phylogeny_lambda_parse_func(pTree ptree, pTreeNode ptnode)
{
	pPhylogenyNode pnode = (pPhylogenyNode)ptnode;
	if (pnode->name) {
		sscanf( pnode->name, "%d", &pnode->taxaid );	
		cafe_param->pcafe->branch_params_cnt++;
	}
	pnode->taxaid--;
}

int __cafe_cmd_lambda_tree(pCafeParam param, char *arg1, char *arg2)
{
	int idx = 1;
	pTree ptree;
	char* plambdastr = NULL;
	param->pcafe->branch_params_cnt = 0;
	if ( arg2 != NULL )
	{
		sscanf( arg1, "%d", &idx );
		plambdastr = arg2;
		ptree = phylogeny_load_from_string(arg2, tree_new, phylogeny_new_empty_node, phylogeny_lambda_parse_func, 0 );
	}
	else
	{
		plambdastr = arg1;
		ptree = phylogeny_load_from_string(arg1, tree_new, phylogeny_new_empty_node, phylogeny_lambda_parse_func, 0 );
	}
	tree_build_node_list(ptree);
	if ( ptree->nlist->size != param->pcafe->super.nlist->size )
	{
		fprintf(stderr, "Lambda has a different topology from the tree\n");
		return -1;
	}
	if (param->pcafe->branch_params_cnt != param->pcafe->super.nlist->size-1) {
		fprintf(stderr,"ERROR(lambda -t): Branch lambda classes not totally specified.\n");
		fprintf(stderr,"%s\n", plambdastr);
		fprintf(stderr,"You have to specify lambda classes for all branches including the internal branches of the tree.\n");
		fprintf(stderr,"There are total %d branches in the tree.\n", param->pcafe->super.nlist->size-1);	// branch_cnt = node_cnt - 1 
		return -1;
	}

	if ( idx == 2 ) 
	{
		if ( tmp_lambda_tree ) phylogeny_free(tmp_lambda_tree);
		tmp_lambda_tree = ptree;
		return 1;
	}
	else 
	{
		if ( param->lambda_tree ) phylogeny_free(param->lambda_tree);
		param->lambda_tree = ptree;
		int l, m, n;
		pArrayList nlist = (pArrayList)param->lambda_tree->nlist;
		memset( param->old_branchlength, 0, sizeof(int) * param->num_branches );	// temporarily use space for old_branchlength 
		for ( l = m = 0 ; l < nlist->size ; l++ )
		{
			int lambda_idx= ((pPhylogenyNode)nlist->array[l])->taxaid;		// lambda tree parameter specification is saved in taxaid
			if ( lambda_idx < 0 ) continue;
			for ( n = 0 ; n < m ; n++ )
			{
				if ( param->old_branchlength[n] == lambda_idx ) break;	// find existing lambda idx
			}
			if ( n == m ) param->old_branchlength[m++] = lambda_idx;	// save new lambda idx
		}
		param->num_lambdas = m;										// number of branch-specific lambdas = m
		if (!param->quiet)
			printf("The number of lambdas is %d\n", m );
	}
	return 0;
}

void set_range_from_family(family_size_range* range, pCafeFamily family)
{
	init_family_size(range, family->max_size);
}

int set_log_file(pCafeParam param, const char *log_file)
{
	if ( param->str_log )
	{
		string_free( param->str_log );
		fclose( param->flog );
		param->str_log = NULL;
	}
	if ( !strcmp(log_file, "stdout" ) )
	{
		param->str_log = NULL;
		param->flog = stdout;
	}
	else
	{
		param->str_log = string_new_with_string(log_file);
		if (  ( param->flog = fopen( param->str_log->buf, "a" ) ) == NULL )
		{
			fprintf(stderr, "ERROR(log): Cannot open log file: %s\n", param->str_log->buf );	
			string_free( param->str_log );
			param->flog = stdout;
			return -1;
		}
	}
	return 0;
}

void __cafe_tree_string_gainloss(pString pstr, pPhylogenyNode ptnode)
{
	int familysize =  ((pCafeNode)ptnode)->familysize;
	if ( ptnode->name ) string_fadd( pstr, "%s", ptnode->name);
	string_fadd(pstr,"_%d", familysize );
	pCafeNode parent = (pCafeNode)ptnode->super.parent;
	if ( parent )
	{
		string_fadd(pstr,"<%d>", familysize - parent->familysize );
	}
}

void __cafe_tree_string_sum_gainloss(pString pstr, pPhylogenyNode ptnode)
{
	int familysize =  ((pCafeNode)ptnode)->familysize;
	if ( ptnode->name ) string_fadd( pstr, "%s", ptnode->name);
	pCafeNode pcnode = (pCafeNode)ptnode;
	string_fadd(pstr,"<%d/%d/%d>", pcnode->viterbi[0], pcnode->viterbi[1], familysize );
}

double __cafe_tree_gainloss_mp_annotation(pString pstr, pTreeNode pnode, pMetapostConfig pmc, va_list ap)
{
	pCafeNode pcnode = (pCafeNode)pnode;
    string_add( pstr, ";\n");
	string_fadd( pstr, "label.urt( btex \\small{%d/%d/%d} ", pcnode->viterbi[0], pcnode->viterbi[1], pcnode->familysize, pnode->id );
	string_fadd( pstr, "etex, p[%d]);\n", pnode->id );
	double last = 0;
	if ( pnode->parent )
	{
		string_fadd( pstr, "xpart mid[%d] = xpart(p[%d]);\n", pnode->id, pnode->id );
		string_fadd( pstr, "ypart mid[%d] = (ypart(p[%d])+ypart(p[%d]))/2;\n", pnode->id, pnode->id, pnode->parent->id );
		string_fadd( pstr, "label.rt( btex $l = %g$ ", ((pPhylogenyNode)pnode)->branchlength );
		string_fadd( pstr, "etex, mid[%d]);\n", pnode->id  );
		string_fadd( pstr, "label.rt( btex $\\lambda=%f$ ", pcnode->birth_death_probabilities.lambda );
		last -= 0.15;
		string_fadd( pstr, "etex, mid[%d] + (0,%fu));\n",  pnode->id, last );
	}
	return last;
}

extern double cafe_tree_mp_remark(pString pstr, pTree ptree, pMetapostConfig pmc, va_list ap1);

int __cafe_cmd_extinct_count_zero(pTree pcafe)
{
	int n;
	int cnt_zero = 0;
	tree_clear_reg(pcafe);
	pArrayList nlist = pcafe->nlist;
	pcafe->root->reg = 1;
	for ( n = 0 ; n < nlist->size ; n+=2 )
	{
		pCafeNode pcnode = (pCafeNode)nlist->array[n];
		if ( pcnode->familysize )
		{
			pTreeNode ptnode = (pTreeNode)pcnode;
			while( ptnode )
			{
				ptnode->reg = 1;
				ptnode = ptnode->parent;
			}
		}
	}
	for ( n = 0 ; n < nlist->size ; n++ )
	{
		pTreeNode pnode = (pTreeNode)nlist->array[n];
		if ( pnode->parent == NULL ) continue;
		if ( pnode->reg == 0 &&  pnode->parent->reg == 1 )
		{
			cnt_zero++;
		}
	}
	return cnt_zero;
}

double __hg_norm_cdf_func(double p, double* args)
{
	return normcdf(p, args[0], args[1]);
}

void __hg_print_sim_extinct(pHistogram** phist_sim_n, pHistogram* phist_sim,  
		                    int r, pHistogram phist_tmp, double* cnt, int num_trials)
{
	int j, t;
	double args[2];
	double alpha;
	for ( j = 0 ; j < phist_sim[r]->nbins ; j++ )
	{
		for ( t = 0 ; t < num_trials ; t++ )
		{
			cnt[t] = histogram_get_count(phist_sim_n[t][r], phist_sim[r]->point[j]);
		}
		histogram_set_by_unit( phist_tmp, cnt, num_trials, 1 );
		args[0] = mean(cnt,num_trials);
		args[1] = sqrt(variance(cnt,num_trials));
		cafe_log(cafe_param, "%g\t%d\t%4.3f\t%g\t%g\t%g ~ %g\n",
				phist_sim[r]->point[j], phist_sim[r]->count[j], ((double)phist_sim[r]->count[j])/phist_sim[r]->nsamples,
				args[0], args[1], phist_tmp->min, phist_tmp->max );		
	}
	memset(cnt, 0, sizeof(double)*num_trials);
	double sum = 0;
	for ( j = 0 ; j < phist_sim[r]->nbins ; j++ )
	{
		double p = phist_sim[r]->point[j];
		if ( p == 0 ) continue;
		for ( t = 0 ; t < num_trials ; t++ )
		{
			double a = p * histogram_get_count(phist_sim_n[t][r], p );
			cnt[t] += a;
			sum += a;
		}
	}
	histogram_set_by_unit( phist_tmp, cnt, num_trials, 1);
	if ( phist_tmp->nbins > 10 )
	{
		histogram_set_by_bin( phist_tmp, cnt, num_trials, 10 );
	}
	args[0] = mean(cnt,num_trials);
	args[1] = sqrt(variance(cnt,num_trials));
	alpha = histogram_check_fitness( phist_tmp, args, __hg_norm_cdf_func );
	cafe_log(cafe_param, "Extinct: %g\t%g\t%g\t%g\t%g ~ %g\n", 
			sum, args[0], args[1], alpha,
			args[0] - 1.96 * args[1], 
			args[0] + 1.96 * args[1]);
	histogram_print( phist_tmp, cafe_param->flog );
}

void cafe_shell_free_errorstruct(pErrorStruct errormodel);


int cafe_shell_read_error_true_measure(const char* errorfile, const char* truefile, int** observed_pairs, int maxFamilySize)
{
	int i=0;
    char buf1[STRING_BUF_SIZE];
    char buf2[STRING_BUF_SIZE];
    
    FILE* fperror = fopen(errorfile,"r");
    if ( fperror == NULL )
    {
        fprintf( stderr, "Cannot open file: %s\n", errorfile );
        return -1;
    }
    if ( fgets(buf1,STRING_BUF_SIZE,fperror) == NULL )
    {
        fclose(fperror);
        fprintf( stderr, "Empty file: %s\n", errorfile );
        return -1;
    }
    FILE* fptruth = fopen(truefile,"r");
    if ( fptruth == NULL )
    {
        fprintf( stderr, "Cannot open file: %s\n", truefile );
        return -1;
    }
    if ( fgets(buf2,STRING_BUF_SIZE,fptruth) == NULL )
    {
        fclose(fptruth);
        fprintf( stderr, "Empty file: %s\n", truefile );
        return -1;
    }
    
    
    // now compare two files and count pairs.
    while(fgets(buf1,STRING_BUF_SIZE,fperror))	
    {
        if ( fgets(buf2,STRING_BUF_SIZE,fptruth) != NULL ) {
            string_pchar_chomp(buf1);
            pArrayList data1 = string_pchar_split( buf1, '\t');
            string_pchar_chomp(buf2);
            pArrayList data2 = string_pchar_split( buf2, '\t');
            if (strcmp((char*)data1->array[1], (char*)data2->array[1])!= 0) {
                fprintf(stderr,"ERROR: the family IDs in each line do not match between the two files\n");      
                return -1;
            }
            // check pairs
            for (i=2; i<data1->size; i++) {
                int size1 = atoi((char*)data1->array[i]);
                int size2 = atoi((char*)data2->array[i]);
                observed_pairs[size1][size2]++;
            }
            arraylist_free(data1,NULL);
            arraylist_free(data2,NULL);
        }
    }
    return 0;
}

// conditional probability of measuring i=familysize when true count is j
int __check_error_model_columnsums(pErrorStruct errormodel)
{
	int i, j = 0;
	int diff = errormodel->todiff;

	for (j = 0; j<diff; j++) {
		// column j
		double columnsum = 0;
		for (i = 0; i <= errormodel->maxfamilysize; i++) {
			columnsum += errormodel->errormatrix[i][j];
		}
		errormodel->errormatrix[0][j] = errormodel->errormatrix[0][j] + (1 - columnsum);
	}

	// all other columns
	for (j = diff; j <= errormodel->maxfamilysize - diff; j++) {
		double columnsum = 0;
		for (i = 0; i <= errormodel->maxfamilysize; i++) {
			columnsum += errormodel->errormatrix[i][j];
		}
		if (abs(1 - columnsum) > 0.00000000000001) {
			for (i = 0; i <= errormodel->maxfamilysize; i++) {
				errormodel->errormatrix[i][j] = errormodel->errormatrix[i][j] / columnsum;
			}
		}
	}

	for (j = errormodel->maxfamilysize - diff + 1; j <= errormodel->maxfamilysize; j++) {
		// column j
		double columnsum = 0;
		for (i = 0; i <= errormodel->maxfamilysize; i++) {
			columnsum += errormodel->errormatrix[i][j];
		}
		errormodel->errormatrix[errormodel->maxfamilysize][j] = errormodel->errormatrix[errormodel->maxfamilysize][j] + (1 - columnsum);

	}
	return 0;
}

pErrorStruct cafe_shell_create_error_matrix_from_estimate(pErrorMeasure errormeasure)
{
	int i = 0;
	int j = 0;

	// allocate new errormodel
	pErrorStruct errormodel = memory_new(1, sizeof(ErrorStruct));

	errormodel->maxfamilysize = errormeasure->maxFamilySize;
	errormodel->fromdiff = -(errormeasure->model_parameter_diff);
	errormodel->todiff = errormeasure->model_parameter_diff;
	errormodel->errorfilename = NULL;
	errormodel->errormatrix = (double**)memory_new_2dim(errormodel->maxfamilysize + 1, errormodel->maxfamilysize + 1, sizeof(double));

	int total_param_num = 0;
	double* total_params = NULL;
	if (errormeasure->b_symmetric) {
		// symmetric
		double sum = 0;
		total_param_num = errormeasure->model_parameter_number + errormeasure->model_parameter_diff + 1;
		total_params = memory_new(total_param_num, sizeof(double));
		total_params[errormeasure->model_parameter_diff] = errormeasure->estimates[0];
		sum = errormeasure->estimates[0];
		for (i = 1; i<errormeasure->model_parameter_number; i++) {
			total_params[errormeasure->model_parameter_diff + i] = errormeasure->estimates[i];
			sum += 2 * errormeasure->estimates[i];
		}
		total_params[total_param_num - 1] = (1 - sum) / (double)((errormeasure->maxFamilySize + 1) - (errormeasure->model_parameter_diff * 2 + 1));
		// now fill left side
		for (i = 0; i<errormeasure->model_parameter_diff; i++) {
			total_params[i] = total_params[abs(total_param_num - 1 - 1 - i)];
		}
	}
	else {
		//asymmetric
		double sum = 0;
		total_param_num = errormeasure->model_parameter_number + 1;
		total_params = memory_new(total_param_num, sizeof(double));
		for (i = 0; i<errormeasure->model_parameter_number; i++) {
			total_params[i] = errormeasure->estimates[i];
			sum += errormeasure->estimates[i];
		}
		total_params[total_param_num - 1] = (1 - sum) / (double)((errormeasure->maxFamilySize + 1) - (errormeasure->model_parameter_diff * 2 + 1));
	}


	// now fill the error matrix column by column
	for (j = 0; j <= errormodel->maxfamilysize; j++) {
		int k = 0;  // k is the index of total_params
		for (i = 0; i<errormodel->fromdiff + j; i++) {
			if (i <= errormodel->maxfamilysize) {
				errormodel->errormatrix[i][j] = total_params[total_param_num - 1]; // marginal error probability epsilon
			}
		}
		for (i = errormodel->fromdiff + j; i <= errormodel->todiff + j; i++) {
			if (i >= 0 && i <= errormodel->maxfamilysize) {
				errormodel->errormatrix[i][j] = total_params[k]; // conditional probability of measuring i+j when true count is j
			}
			k++;
		}
		for (i = errormodel->todiff + j + 1; i <= errormodel->maxfamilysize; i++) {
			if (i >= 0) {
				errormodel->errormatrix[i][j] = total_params[total_param_num - 1]; // marginal error probability epsilon
			}
		}
	}

	// now make sure that columns of the error matrix sums to one.
	__check_error_model_columnsums(errormodel);

	/*    if (cafe_param->pfamily->errors == NULL) {
	cafe_param->pfamily->errors = arraylist_new(cafe_param->pfamily->num_species);
	}
	arraylist_add(cafe_param->pfamily->errors, errormodel);*/
	return errormodel;
}

 



double __loglikelihood_pairs_from_true_measure(double* parameters, void* args)
{
	int i, j;
    
    pErrorMeasure errormeasure = (pErrorMeasure)args;
    
    double marginal_error_probability_epsilon = 0;   
    if (errormeasure->b_symmetric) { 
        // symmetric
        double sum = parameters[0];
        for(i=1; i<errormeasure->model_parameter_number; i++) {
            sum += 2*parameters[i];
        } 
        marginal_error_probability_epsilon = (1-sum)/(double)((errormeasure->maxFamilySize+1)-(errormeasure->model_parameter_diff*2+1));
    }   
    else {  
        //asymmetric
        double sum = 0;
        for(i=0; i<errormeasure->model_parameter_number; i++) {
            sum += parameters[i];
        } 
        marginal_error_probability_epsilon = (1-sum)/(double)((errormeasure->maxFamilySize+1)-(errormeasure->model_parameter_diff*2+1));
    }
  
    
	double score = 0;
	int skip = 0;
	for ( i = 0 ; i < errormeasure->model_parameter_number ; i++ )
	{
		if ( ( parameters[i] < 0 ) || ( marginal_error_probability_epsilon < 0 ) || (marginal_error_probability_epsilon > parameters[i]) )
		{ 
			skip  = 1;
			score = log(0);
			break;
		}
	}
	if ( !skip && errormeasure->b_peakzero ) {
        double previous_param = 0;
        if (errormeasure->b_symmetric) {
            previous_param = parameters[0]; 
            for (i = 1; i<errormeasure->model_parameter_number; i++) {
                if (previous_param < parameters[i]) {
                    skip  = 1;
                    score = log(0);
                    break;                
                }
                previous_param = parameters[i];
            }
        }
        else {
            previous_param = parameters[errormeasure->model_parameter_diff]; 
            for (i=1; i<=errormeasure->model_parameter_diff; i++) {
                if (previous_param < parameters[errormeasure->model_parameter_diff-i]) {
                    skip  = 1;
                    score = log(0);
                    break;                
                }
                previous_param = parameters[errormeasure->model_parameter_diff-i];
            }
            previous_param = parameters[errormeasure->model_parameter_diff]; 
            for (i=1; i<=errormeasure->model_parameter_diff; i++) {
                if (previous_param < parameters[errormeasure->model_parameter_diff+i]) {
                    skip  = 1;
                    score = log(0);
                    break;                
                }
                previous_param = parameters[errormeasure->model_parameter_diff+i];
            }
        }
    }
	if ( !skip )
	{
        errormeasure->estimates = parameters;
        pErrorStruct errormodel = cafe_shell_create_error_matrix_from_estimate(errormeasure);
        
        double** discord_prob_model = (double**)memory_new_2dim(errormeasure->maxFamilySize+1, errormeasure->maxFamilySize+1, sizeof(double));
        for (i=0; i<=errormeasure->maxFamilySize; i++) {
            for (j=0; j<=errormeasure->maxFamilySize; j++) {
                // find discordance probability based on parameters
                double pi_i_j = errormodel->errormatrix[i][j];
                discord_prob_model[i][j] = errormeasure->sizeDist[j]*pi_i_j;
            }
        }
        for (i=0; i<=errormeasure->maxFamilySize; i++) {
            for (j=0; j<=errormeasure->maxFamilySize; j++) {
                // add to the log likelihood
                double term = errormeasure->pairs[i][j]? errormeasure->pairs[i][j] * log(discord_prob_model[i][j]) : 0;
                score += term;
                if (isnan(score) || isinf(-score)) {
                    cafe_log(cafe_param,"Score: %f\n", score);
                }
            }
        }
        double prob00 = errormodel->errormatrix[0][0]*errormeasure->sizeDist[0];
        score -= log(1-prob00);

        memory_free_2dim((void**)discord_prob_model, errormeasure->maxFamilySize+1, errormeasure->maxFamilySize+1, NULL);
        cafe_shell_free_errorstruct(errormodel);
        
    }
    
    char buf[STRING_STEP_SIZE];
	buf[0] = '\0';
	string_pchar_join_double(buf,",", errormeasure->model_parameter_number, parameters );
	cafe_log(cafe_param,"\tparameters : %s & Score: %f\n", buf, score);
    return -score;
}





void cafe_shell_free_errorstruct(pErrorStruct errormodel)
{
    if (errormodel->errorfilename) {
        memory_free(errormodel->errorfilename);
        errormodel->errorfilename = NULL;
    }
    if (errormodel->errormatrix) {
        memory_free_2dim((void**)errormodel->errormatrix, errormodel->maxfamilysize+1, errormodel->maxfamilysize+1, NULL);  
        errormodel->errormatrix = NULL;
    }
}

