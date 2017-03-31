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

extern pBirthDeathCacheArray probability_cache;

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

void cafe_shell_prompt(char* prompt, char* format, ... )
{
	va_list ap;
	va_start(ap, format);
	printf("%s ", prompt);
	if (vscanf( format, ap ) == EOF)
		fprintf(stderr, "Read failure\n");
	va_end(ap);
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

int cafe_shell_set_familysize()
{
	int i;
	int max = 0;
	char buf[STRING_STEP_SIZE];

	STDERR_IF( cafe_param->pcafe == NULL, "You did not specify tree: command 'tree'\n" );

	pArrayList nlist = cafe_param->pcafe->super.nlist;
	for ( i = 0; i < nlist->size ; i+=2 )
	{
		pCafeNode pnode = (pCafeNode)nlist->array[i];
		sprintf(buf, "%s: ", pnode->super.name );
		int size = -1;
		cafe_shell_prompt( buf , "%d", &size  );
		if ( size < 0 )
		{ 
			fprintf( stderr, "ERROR: You put wrong data, you must enter an integer greater than or equal to 0\n");
			cafe_shell_prompt( "Retry? [Y|N] ", "%s", buf);
			if ( buf[0] != 'Y' && buf[0] != 'y' ) return -1; 
			i -= 2;
		}
		else
		{
			pnode->familysize = size;
			if ( size > max ) max = size;
		}
	}
	return max;
}

int cafe_shell_set_branchlength()
{
	int i;
	char buf[STRING_STEP_SIZE];

	pArrayList nlist = cafe_param->pcafe->super.nlist;
	for ( i = 0; i < nlist->size ; i++ )
	{
		pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
		if ( tree_is_root( (pTree)cafe_param->pcafe, (pTreeNode)pnode) ) continue;
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
	if (probability_cache) cafe_tree_set_birthdeath(cafe_param->pcafe, probability_cache);
	return 0;
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

int __cafe_cmd_lambda_tree(pArgument parg)
{
	int idx = 1;
	pTree ptree;
	char* plambdastr = NULL;
	cafe_param->pcafe->branch_params_cnt = 0;
	if ( parg->argc == 2 )
	{
		sscanf( parg->argv[0], "%d", &idx );
		plambdastr = parg->argv[1];
		ptree = phylogeny_load_from_string(parg->argv[1], tree_new, phylogeny_new_empty_node, phylogeny_lambda_parse_func, 0 );
	}
	else
	{
		plambdastr = parg->argv[0];
		ptree = phylogeny_load_from_string(parg->argv[0], tree_new, phylogeny_new_empty_node, phylogeny_lambda_parse_func, 0 );
	}
	tree_build_node_list(ptree);
	if ( ptree->nlist->size != cafe_param->pcafe->super.nlist->size )
	{
		fprintf(stderr, "Lambda has a different topology from the tree\n");
		return -1;
	}
	if (cafe_param->pcafe->branch_params_cnt != cafe_param->pcafe->super.nlist->size-1) {
		fprintf(stderr,"ERROR(lambda -t): Branch lambda classes not totally specified.\n");
		fprintf(stderr,"%s\n", plambdastr);
		fprintf(stderr,"You have to specify lambda classes for all branches including the internal branches of the tree.\n");
		fprintf(stderr,"There are total %d branches in the tree.\n", cafe_param->pcafe->super.nlist->size-1);	// branch_cnt = node_cnt - 1 
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
		if ( cafe_param->lambda_tree ) phylogeny_free(cafe_param->lambda_tree);
		cafe_param->lambda_tree = ptree;
		int l, m, n;
		pArrayList nlist = (pArrayList)cafe_param->lambda_tree->nlist;
		memset( cafe_param->old_branchlength, 0, sizeof(int) * cafe_param->num_branches );	// temporarily use space for old_branchlength 
		for ( l = m = 0 ; l < nlist->size ; l++ )
		{
			int lambda_idx= ((pPhylogenyNode)nlist->array[l])->taxaid;		// lambda tree parameter specification is saved in taxaid
			if ( lambda_idx < 0 ) continue;
			for ( n = 0 ; n < m ; n++ )
			{
				if ( cafe_param->old_branchlength[n] == lambda_idx ) break;	// find existing lambda idx
			}
			if ( n == m ) cafe_param->old_branchlength[m++] = lambda_idx;	// save new lambda idx
		}
		cafe_param->num_lambdas = m;										// number of branch-specific lambdas = m
		if (!cafe_param->quiet)
			printf("The number of lambdas is %d\n", m );
	}
	return 0;
}








void __cafe_cmd_viterbi_family_print(int idx)
{
	pCafeTree pcafe = cafe_param->pcafe;
	cafe_family_set_size_with_family_forced(cafe_param->pfamily,idx,pcafe);
	compute_tree_likelihoods(pcafe);
	int ridx =  __maxidx(((pCafeNode)pcafe->super.root)->likelihoods,pcafe->rfsize) + pcafe->rootfamilysizes[0];
	double mlh =  __max( ((pCafeNode)pcafe->super.root)->likelihoods,pcafe->rfsize);
	//compute_tree_likelihoods(pcafe);
	cafe_tree_viterbi(pcafe);
	pString pstr = cafe_tree_string(pcafe);
	printf("%g(%d)\t%s\n", mlh , ridx,  pstr->buf );
	string_free(pstr);
}





double _cafe_cross_validate_by_family(const char* queryfile, const char* truthfile, const char* errortype) 
{
	int i, j;
	double MSE = 0;
	double MAE = 0;
	double SSE = 0;
	double SAE = 0;
	cafe_family_read_query_family(cafe_param, queryfile);
	if ( cafe_param->cv_test_count_list == NULL ) return -1;
	
	// read in validation data
	pCafeFamily truthfamily = cafe_family_new( truthfile, 1 );
	if ( truthfamily == NULL ) {
		fprintf(stderr, "failed to read in true values %s\n", truthfile);
		return -1;
	}
	
	// now compare reconstructed count to true count	
	pCafeTree pcafe = cafe_param->pcafe;
	pCafeTree truthtree = cafe_tree_copy(pcafe);
	// set parameters
	if ( truthtree )
	{
		cafe_family_set_species_index(truthfamily, truthtree);
	}

	reset_birthdeath_cache(cafe_param->pcafe, cafe_param->parameterized_k_value, &cafe_param->family_size);
	
	for(i=0; i< cafe_param->cv_test_count_list->size; i++) 
	{
		int* testcnt = (int*)cafe_param->cv_test_count_list->array[i];
		cafe_family_set_size(truthfamily, i, truthtree);
		cafe_family_set_size_by_species(cafe_param->cv_test_species_list->array[i], *testcnt, pcafe);
		if (cafe_param->posterior) {
			cafe_tree_viterbi_posterior(pcafe, cafe_param);
		}
		else {
			cafe_tree_viterbi(pcafe);
		}
		// leaf nodes SSE
		SSE = 0;
		SAE = 0;
		int nodecnt = 0;
		for (j=0; j<pcafe->super.nlist->size; j=j+2) {
			int error = ((pCafeNode)truthtree->super.nlist->array[j])->familysize-((pCafeNode)pcafe->super.nlist->array[j])->familysize;
			SSE += pow(error, 2);
			SAE += abs(error);
			nodecnt++;
		}
		MSE += SSE/nodecnt;
		MSE += SAE/nodecnt;
	}
	cafe_free_birthdeath_cache(pcafe);

	MSE = MSE/cafe_param->cv_test_count_list->size;
	MAE = MAE/cafe_param->cv_test_count_list->size;
	cafe_log( cafe_param, "MSE %f\n", MSE );
	cafe_log( cafe_param, "MAE %f\n", MSE );

	double returnerror = -1;
	if (strncmp(errortype, "MSE", 3)==0) {
		returnerror = MSE;
	}
	else if (strncmp(errortype, "MAE", 3)==0) {
		returnerror = MAE;
	}
	return returnerror;
}




double _cafe_cross_validate_by_species(const char* validatefile, const char* errortype) 
{
	int i, j;
	cafe_family_read_validate_species( cafe_param, validatefile );
	if ( cafe_param->cv_test_count_list == NULL ) return -1;
	// now compare reconstructed count to true count	
	pCafeTree pcafe = cafe_param->pcafe;

	reset_birthdeath_cache(cafe_param->pcafe, cafe_param->parameterized_k_value, &cafe_param->family_size);
	pArrayList estimate_size = arraylist_new(cafe_param->cv_test_count_list->size);
	for(i=0; i< cafe_param->pfamily->flist->size; i++) 
	{
		cafe_family_set_size(cafe_param->pfamily,i, pcafe);
		if (cafe_param->posterior) {
			cafe_tree_viterbi_posterior(pcafe, cafe_param);
		}
		else {
			cafe_tree_viterbi(pcafe);
		}
		for (j=0; j<pcafe->super.nlist->size; j++) {
			char* nodename = ((pPhylogenyNode)pcafe->super.nlist->array[j])->name;
			if (nodename && (strcmp(nodename, cafe_param->cv_species_name)==0)) {
				int* pFamilysize = memory_new(1, sizeof(int));
				*pFamilysize = ((pCafeNode)pcafe->super.nlist->array[j])->familysize;
				arraylist_add(estimate_size, (void*)pFamilysize);
			}
			
		}
	}
	cafe_free_birthdeath_cache(pcafe);
	double MSE = 0;
	double MAE = 0;
	STDERR_IF(cafe_param->cv_test_count_list->size != cafe_param->pfamily->flist->size, "list size don't match\n");
	for(i=0; i<cafe_param->cv_test_count_list->size; i++) {
		int error = (*((int*)cafe_param->cv_test_count_list->array[i]) - *((int*)estimate_size->array[i]));
		MSE += pow(error, 2);
		MAE += abs(error);
	}
	MSE = MSE/(cafe_param->cv_test_count_list->size);
	MAE = MAE/(cafe_param->cv_test_count_list->size);
	cafe_log( cafe_param, "MSE %f\n", MSE );
	cafe_log( cafe_param, "MAE %f\n", MAE );
	
	arraylist_free(estimate_size, free);
	
	double returnerror = -1;
	if (strncmp(errortype, "MSE", 3)==0) {
		returnerror = MSE;
	}
	else if (strncmp(errortype, "MAE", 3)==0) {
		returnerror = MAE;
	}
	return returnerror;
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


int cafe_shell_read_freq_from_measures(const char* file1, const char* file2, int* sizeFreq)
{
    int i=0;
    char buf1[STRING_BUF_SIZE];
    char buf2[STRING_BUF_SIZE];
    
    FILE* fpfile1 = fopen(file1,"r");
    if ( fpfile1 == NULL )
    {
        fprintf( stderr, "Cannot open file: %s\n", file1 );
        return -1;
    }
    if ( fgets(buf1,STRING_BUF_SIZE,fpfile1) == NULL )
    {
        fclose(fpfile1);
        fprintf( stderr, "Empty file: %s\n", file1 );
        return -1;
    }
    FILE* fpfile2 = NULL;
    if (file2) {
        fpfile2 = fopen(file2,"r");
        if ( fpfile2 == NULL )
        {
            fprintf( stderr, "Cannot open file: %s\n", file2 );
            return -1;
        }
        if ( fgets(buf2,STRING_BUF_SIZE,fpfile2) == NULL )
        {
            fclose(fpfile2);
            fprintf( stderr, "Empty file: %s\n", file2 );
            return -1;
        }
    }
    
    string_pchar_chomp(buf1);
    pArrayList data1 = string_pchar_split( buf1, '\t');
    arraylist_free(data1,NULL);        
    
    // count frequency of family sizes
    int line1 = 0;
    int maxFamilySize = cafe_param->family_size.max;
    int data1colnum = 0;
    while(fgets(buf1,STRING_BUF_SIZE,fpfile1))	
    {
        string_pchar_chomp(buf1);
        pArrayList data1 = string_pchar_split( buf1, '\t');
        for (i=2; i<data1->size; i++) {
            int size1 = atoi((char*)data1->array[i]);
            sizeFreq[size1]++;
            if (size1 > maxFamilySize) {
                maxFamilySize = size1;
            }
        }
        data1colnum = data1->size;
        arraylist_free(data1,NULL);
        line1++;
    }
    if (fpfile2) 
    {
        int line2 = 0;
        string_pchar_chomp(buf2);
        pArrayList data2 = string_pchar_split( buf2, '\t');
        if (data1colnum != data2->size) {
            fprintf(stderr,"file: the number of columns do not match between the two files\n");
            return -1;
        }
        arraylist_free(data2, NULL);

        while(fgets(buf2,STRING_BUF_SIZE,fpfile2))	
        {
            string_pchar_chomp(buf2);
            pArrayList data2 = string_pchar_split( buf2, '\t');
            for (i=2; i<data2->size; i++) {
                int size2 = atoi((char*)data2->array[i]);
                sizeFreq[size2]++;
                if (size2 > maxFamilySize) {
                    maxFamilySize = size2;
                }
            }
            arraylist_free(data2,NULL);
            line2++;
        }    
        if (line1 != line2) {
            fprintf(stderr,"ERROR: the number of lines do not match between the two files\n");
            return -1;
        }
    }
    return maxFamilySize;
}




int cafe_shell_read_error_double_measure(const char* error1, const char* error2, int** observed_pairs, int maxFamilySize)
{
	int i=0;
    int j=0;
    char buf1[STRING_BUF_SIZE];
    char buf2[STRING_BUF_SIZE];
    
    FILE* fperror1 = fopen(error1,"r");
    if ( fperror1 == NULL )
    {
        fprintf( stderr, "Cannot open file: %s\n", error1 );
        return -1;
    }
    if ( fgets(buf1,STRING_BUF_SIZE,fperror1) == NULL )
    {
        fclose(fperror1);
        fprintf( stderr, "Empty file: %s\n", error1 );
        return -1;
    }
    FILE* fperror2 = fopen(error2,"r");
    if ( fperror2 == NULL )
    {
        fprintf( stderr, "Cannot open file: %s\n", error2 );
        return -1;
    }
    if ( fgets(buf2,STRING_BUF_SIZE,fperror2) == NULL )
    {
        fclose(fperror2);
        fprintf( stderr, "Empty file: %s\n", error2 );
        return -1;
    }
    
  
    // now compare two files and count pairs.
    while(fgets(buf1,STRING_BUF_SIZE,fperror1))	
    {
        if ( fgets(buf2,STRING_BUF_SIZE,fperror2) != NULL ) {
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
    
    // now make triangle matrix by merging i,j and j,i
    for (i=0; i<=maxFamilySize; i++) {
        for (j=0; j<i; j++) {
            observed_pairs[j][i] += observed_pairs[i][j];
            observed_pairs[i][j] = 0;
        }
    }
    
    
    return 0;
}

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

double __loglikelihood_pairs_from_double_measure(double* parameters, void* args)
{
	int i, j, k;
    
    pErrorMeasure errormeasure = (pErrorMeasure)args;
    double marginal_error_probability_epsilon = -1;   
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
            for (j=i; j<=errormeasure->maxFamilySize; j++) {
                for (k=0; k<=errormeasure->maxFamilySize; k++) {
                    double pi_i_k = errormodel->errormatrix[i][k];
                    double pi_j_k = errormodel->errormatrix[j][k];
                    if (i==j) {
                        discord_prob_model[i][j] += errormeasure->sizeDist[k]*pi_i_k*pi_j_k;
                    }
                    else {
                        discord_prob_model[i][j] += 2*errormeasure->sizeDist[k]*pi_i_k*pi_j_k;                        
                    }
                }
            }
        }
        for (i=0; i<=errormeasure->maxFamilySize; i++) {
            for (j=i; j<=errormeasure->maxFamilySize; j++) {
                // add to the log likelihood
                double term = errormeasure->pairs[i][j]? errormeasure->pairs[i][j] * log(discord_prob_model[i][j]) : 0;
                score += term;
                if (isnan(score) || isinf(-score) || !isfinite(score)) {
                    cafe_log(cafe_param,"Score: %f\n", score);
                    break;
                }
            }
        }
        double prob00 = 0;
        for (k=0; k<=errormeasure->maxFamilySize; k++) {
            double pi_i_k = errormodel->errormatrix[0][k];
            double pi_j_k = errormodel->errormatrix[0][k];
            prob00 += errormeasure->sizeDist[k]*pi_i_k*pi_j_k;
        }
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




pErrorMeasure cafe_shell_estimate_error_double_measure(const char* error1, const char* error2, int b_symmetric, int max_diff, int b_peakzero)
{
    int i;
    pCafeParam param = cafe_param;
    
    int* sizeFreq = memory_new(10000, sizeof(int)); 
    int maxFamilySize = cafe_shell_read_freq_from_measures(error1, error2, sizeFreq);
    if (maxFamilySize < 0 ) {
        fprintf(stderr,"ERROR: failed to read freqeuncy from measurement files\n");              
    }
    // get size probability distribution
    int sizeTotal = 0;
    for (i = 0; i<= maxFamilySize; i++) {
        sizeTotal += sizeFreq[i]+1;
        if (sizeTotal < 0) {
            fprintf(stderr,"ERROR: total freqeuncy is less than zero\n");                          
        }
    }
    double* sizeDist = (double*)memory_new(maxFamilySize+1, sizeof(double));
    for (i = 0; i<=maxFamilySize; i++) {
        sizeDist[i] = (sizeFreq[i]+1)/(double)sizeTotal;
        if (sizeDist[i] < 0) {
            fprintf(stderr,"ERROR: freqeuncy is less than zero\n");                          
        }
    }
    
    
    int** observed_pairs = (int**)memory_new_2dim(maxFamilySize+1, maxFamilySize+1, sizeof(int)); // need space for zero
    int retval = cafe_shell_read_error_double_measure(error1, error2, observed_pairs, maxFamilySize);
    if (retval < 0) {
        fprintf(stderr,"ERROR: failed to count pairs from measurement files\n");      
    }
    
    // set up parameters for ML
    pErrorMeasure error = memory_new(1, sizeof(ErrorMeasure));
    error->sizeDist = sizeDist;
    error->maxFamilySize = maxFamilySize;
    error->pairs = observed_pairs;
    error->b_symmetric = b_symmetric;
    error->b_peakzero = b_peakzero;
    if (b_symmetric) {
        // symmetric model (diff == number)
        error->model_parameter_diff = max_diff;
        error->model_parameter_number = max_diff+1;  
    }
    else {
        // asymmetric model (diff*2 == number)
        error->model_parameter_diff = max_diff;
        error->model_parameter_number = 2*max_diff+1;
    }
    
        
    // now estimate the misclassification rate 
    int max_runs = 100;
    int converged = 0;
	int runs = 0;
    double minscore = DBL_MAX;
    double* parameters = memory_new(error->model_parameter_number, sizeof(double));
	double* bestrun_parameters = memory_new(error->model_parameter_number, sizeof(double));
    
    do {
        pFMinSearch pfm;
        double* sorted_params = memory_new_with_init(error->model_parameter_number, sizeof(double), parameters);
        for (i=0; i<error->model_parameter_number; i++) {
            sorted_params[i] = unifrnd()/(double)error->model_parameter_number;
        }
        qsort (sorted_params, error->model_parameter_number, sizeof(double), comp_double);
        if (error->b_symmetric) {
            int j=0;
            for (i=error->model_parameter_number-1; i>=0; i--) {
                parameters[j++] = sorted_params[i];
            }
        }
        else {
            int j=error->model_parameter_number-1;
            parameters[error->model_parameter_diff] = sorted_params[j--];
            for (i=1; i<=error->model_parameter_diff; i++) {
                parameters[error->model_parameter_diff-i] = sorted_params[j--];
                parameters[error->model_parameter_diff+i] = sorted_params[j--];
            }
        }
        pfm = fminsearch_new_with_eq(__loglikelihood_pairs_from_double_measure, error->model_parameter_number, error);
        pfm->tolx = 1e-9;
        pfm->tolf = 1e-9;
        fminsearch_min(pfm, parameters);
        double *re = fminsearch_get_minX(pfm);
        for ( i = 0 ; i < error->model_parameter_number; i++ ) parameters[i] = re[i];
        cafe_log(param, "\n");
        cafe_log(param,"Misclassification Matrix Search Result: %d\n", pfm->iters );
        cafe_log(param, "Score: %f\n", *pfm->fv);
        
        if (runs > 0) {
			if (!isnan(*pfm->fv) && !isinf(*pfm->fv) && abs(minscore - (*pfm->fv)) < pfm->tolf) {
				converged = 1;
			}
		}
        if (pfm->iters < pfm->maxiters) {
            if ( *pfm->fv < minscore) {
                minscore = *pfm->fv;
                memcpy(bestrun_parameters, parameters, (error->model_parameter_number)*sizeof(double));
            }
            runs++;
        }
/*        else {
            cafe_log(param,"what went wrong?\n");
            fminsearch_min(pfm, parameters);
        }*/
		fminsearch_free(pfm);
	} while (!converged && runs<max_runs); 

		if (converged) {
			cafe_log(param,"score converged in %d runs.\n", runs);
		}
		else {
			cafe_log(param,"score failed to converge in %d runs.\n", max_runs);
			cafe_log(param,"best score: %f\n", minscore);            
		}
	memory_free(parameters);      
    error->estimates = bestrun_parameters;
    
    //memory_free(error);           // we are going to return these values
    memory_free_2dim((void**)observed_pairs, maxFamilySize+1, maxFamilySize+1, NULL);
    memory_free(sizeFreq);
    return error; 
}



pErrorMeasure cafe_shell_estimate_error_true_measure(const char* errorfile, const char* truefile, int b_symmetric, int max_diff, int b_peakzero)
{
    int i;
    pCafeParam param = cafe_param;
    
    int* sizeFreq = memory_new(10000, sizeof(int)); 
    int maxFamilySize = cafe_shell_read_freq_from_measures(truefile, errorfile, sizeFreq);
    if (maxFamilySize < 0 ) {
        fprintf(stderr,"ERROR: failed to read freqeuncy from measurement files\n");              
    }
    // get size probability distribution
    int sizeTotal = 0;
    for (i = 0; i<= maxFamilySize; i++) {
        sizeTotal += sizeFreq[i]+1;
    }
    double* sizeDist = (double*)memory_new(maxFamilySize+1, sizeof(double));
    for (i = 0; i<=maxFamilySize; i++) {
        sizeDist[i] = (sizeFreq[i]+1)/(double)sizeTotal;
    }
    
    
    int** observed_pairs = (int**)memory_new_2dim(maxFamilySize+1, maxFamilySize+1, sizeof(int)); // need space for zero
    int retval = cafe_shell_read_error_true_measure(errorfile, truefile, observed_pairs, maxFamilySize);
    if (retval < 0) {
        fprintf(stderr,"ERROR: failed to count pairs from measurement files\n");      
    }
    
    // set up parameters for ML
    pErrorMeasure error = memory_new(1, sizeof(ErrorMeasure));
    error->sizeDist = sizeDist;
    error->maxFamilySize = maxFamilySize;
    error->pairs = observed_pairs;
    error->b_symmetric = b_symmetric;
    error->b_peakzero = b_peakzero;
    if (b_symmetric) {
        // symmetric model (diff == number)
        error->model_parameter_diff = max_diff;
        error->model_parameter_number = max_diff+1;  
    }
    else {
        // asymmetric model (diff*2 == number)
        error->model_parameter_diff = max_diff;
        error->model_parameter_number = 2*max_diff+1;
    }
    
    // now estimate the misclassification rate 
    int max_runs = 100;
    int converged = 0;
	int runs = 0;
    double minscore = DBL_MAX;
    double* parameters = memory_new(error->model_parameter_number, sizeof(double));
	double* bestrun_parameters = memory_new(error->model_parameter_number, sizeof(double));
    
    do {
        pFMinSearch pfm;
        double* sorted_params = memory_new_with_init(error->model_parameter_number, sizeof(double), parameters);
        for (i=0; i<error->model_parameter_number; i++) {
            sorted_params[i] = unifrnd()/(double)error->model_parameter_number;
        }
        qsort (sorted_params, error->model_parameter_number, sizeof(double), comp_double);
        if (error->b_symmetric) {
            int j=0;
            for (i=error->model_parameter_number-1; i>=0; i--) {
                parameters[j++] = sorted_params[i];
            }
        }
        else {
            int j=error->model_parameter_number-1;
            parameters[error->model_parameter_diff] = sorted_params[j--];
            for (i=1; i<=error->model_parameter_diff; i++) {
                parameters[error->model_parameter_diff-i] = sorted_params[j--];
                parameters[error->model_parameter_diff+i] = sorted_params[j--];
            }
        }
        
        pfm = fminsearch_new_with_eq(__loglikelihood_pairs_from_true_measure, error->model_parameter_number, error);
        pfm->tolx = 1e-9;
        pfm->tolf = 1e-9;
        fminsearch_min(pfm, parameters);
        double *re = fminsearch_get_minX(pfm);
        for ( i = 0 ; i < error->model_parameter_number; i++ ) parameters[i] = re[i];
        cafe_log(param, "\n");
        cafe_log(param,"Misclassification Matrix Search Result: %d\n", pfm->iters );
        cafe_log(param, "Score: %f\n", *pfm->fv);
        
        if (runs > 0) {
			if (!isnan(*pfm->fv) && !isinf(*pfm->fv) && abs(minscore - (*pfm->fv)) < pfm->tolf) {
				converged = 1;
			}
		}
        if (pfm->iters < pfm->maxiters) {
            if ( *pfm->fv < minscore) {
                minscore = *pfm->fv;
                memcpy(bestrun_parameters, parameters, (error->model_parameter_number)*sizeof(double));
            }
            runs++;
        }
		fminsearch_free(pfm);
	} while (!converged && runs<max_runs); 
    
    if (converged) {
        cafe_log(param,"score converged in %d runs.\n", runs);
    }
    else {
        cafe_log(param,"score failed to converge in %d runs.\n", max_runs);
        cafe_log(param,"best score: %f\n", minscore);            
    }
	memory_free(parameters);      
    error->estimates = bestrun_parameters;
    
    //memory_free(error);           // we are going to return these values
    memory_free_2dim((void**)observed_pairs, maxFamilySize+1, maxFamilySize+1, NULL);
    memory_free(sizeFreq);
    return error; 
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

