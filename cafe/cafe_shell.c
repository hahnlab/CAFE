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
#include<io.h>
#include "cafe.h"
#include "cafe_shell.h"
#include "viterbi.h"

extern int cafe_shell_dispatch_commandf(char* format, ...);

/**
* \brief Holds the global program state that user commands act on.
*
*/
pCafeParam cafe_param;

pTree tmp_lambda_tree;
pArrayList cafe_pCD;

#ifndef STDERR_IF
	#define STDERR_IF(a,b)	if ( a ) { fprintf(stderr,b); return -1; }
#endif

/**
* \brief Holds the list of commands that are available in Cafe.
* 
* Each element consists of a command and the function that is 
* called to handle that command. Terminated by a NULL, NULL pair
* Functions include #cafe_cmd_lambda, #cafe_cmd_family, #cafe_cmd_tree,
* etc.
*/
CafeShellCommand cafe_cmd[]  =
{
	{ "branchlength", cafe_cmd_branchlength },
	{ "extinct", cafe_cmd_extinct },
	{ "family", cafe_cmd_family },
	{ "genfamily", cafe_cmd_generate_random_family },
	{ "info", cafe_cmd_print_param },
	{ "lambdamu", cafe_cmd_lambda_mu },
	{ "lhtest", cafe_cmd_lh_test },
	{ "load", cafe_cmd_load },
	{ "log", cafe_cmd_log },
	{ "pvalue", cafe_cmd_pvalue },
	{ "retrieve", cafe_cmd_retrieve },
	{ "rootdist", cafe_cmd_root_dist}, 
    { "errormodel", cafe_cmd_error_model},  
    { "noerrormodel", cafe_cmd_no_error_model},  
    { "simerror", cafe_cmd_simerror},
    { "esterror", cafe_cmd_esterror},
	{ "cvspecies", cafe_cmd_crossvalidation_by_species},
	{ "cvfamily", cafe_cmd_crossvalidation_by_family},
	{ "accuracy", cafe_cmd_reconstruction_accuracy},
	{ "save", cafe_cmd_save},
	{ "score", cafe_cmd_score },
	{ "simextinct", cafe_cmd_sim_extinct },
	{ "tree", cafe_cmd_tree },
	{ "version", cafe_cmd_version },
	{ "viterbi", cafe_cmd_viterbi},
	{ NULL, NULL }
};

int __cafe_cmd_log(int argc ,char* argv[] );


pArrayList cafe_shell_build_argument(int argc, char* argv[])
{
	int i, j;
	pArrayList pal = arraylist_new(20);
	for ( i = 1 ; i < argc ; i++ )
	{
		if ( argv[i][0] == '-' && !isdigit(argv[i][1]) )
		{
			pArgument parg = (pArgument)memory_new(1,sizeof(Argument));
			parg->argc = 0;
			parg->opt = argv[i];
			for ( j = i+1; j < argc; j++ )
			{
				if ( argv[j][0] == '-' && !isdigit(argv[j][1]) ) break;
				parg->argc++;
			}
			parg->argv = parg->argc ? &argv[i+1]  : NULL;
			arraylist_add( pal, parg );
			i = j - 1;
		}
	}
	return pal;
}


pArgument cafe_shell_get_argument(char* opt, pArrayList pal)
{
	int i;
	for ( i = 0 ; i < pal->size; i++ )
	{
		pArgument parg = (pArgument)pal->array[i];
		if ( !strcasecmp(parg->opt, opt) )
		{
			return parg;
		}
	}
	return NULL;
}


void cafe_shell_set_lambda(pCafeParam param, double* lambda);
void cafe_shell_set_lambda_mu(pCafeParam param, double* parameters);

void viterbi_parameters_init(viterbi_parameters *viterbi, int nnodes, int nrows)
{
	viterbi->num_nodes = nnodes;
	viterbi->num_rows = nrows;
	viterbi->viterbiPvalues = (double**)memory_new_2dim(nnodes, nrows, sizeof(double));
	viterbi->expandRemainDecrease = (int**)memory_new_2dim(3, nnodes, sizeof(int));
	viterbi->viterbiNodeFamilysizes = (int**)memory_new_2dim(nnodes, nrows, sizeof(int));
	viterbi->maximumPvalues = (double*)memory_new(nrows, sizeof(double));
	viterbi->averageExpansion = (double*)memory_new(nnodes, sizeof(double));
}

void viterbi_parameters_clear(viterbi_parameters* viterbi, int nnodes)
{
//	viterbi_parameters* viterbi = &param->viterbi;
	if ( viterbi->viterbiPvalues )
	{
		int num = (nnodes - 1 )/2;
		memory_free_2dim((void**)viterbi->viterbiPvalues,num,0,NULL);
		memory_free_2dim((void**)viterbi->expandRemainDecrease,3,0,NULL);
		memory_free_2dim((void**)viterbi->viterbiNodeFamilysizes,num, 0, NULL);
		memory_free(viterbi->averageExpansion);
		viterbi->averageExpansion = NULL;
		if ( viterbi->maximumPvalues )
		{
			memory_free(viterbi->maximumPvalues);
			viterbi->maximumPvalues = NULL;
		}
	}
	if ( viterbi->cutPvalues )
	{
		memory_free_2dim((void**)viterbi->cutPvalues,nnodes,0,NULL);
	}
	viterbi->viterbiPvalues = NULL;
	viterbi->expandRemainDecrease = NULL;
	viterbi->viterbiNodeFamilysizes = NULL;
	viterbi->cutPvalues = NULL;
	viterbi->maximumPvalues = NULL;
}

void viterbi_set_max_pvalue(viterbi_parameters* viterbi, int index, double val)
{
	assert(index < viterbi->num_rows);
	viterbi->maximumPvalues[index] = val;
}


void cafe_shell_clear_param(pCafeParam param, int btree_skip)
{
	if ( param->str_fdata ) 
	{
		string_free( param->str_fdata );
		param->str_fdata = NULL;
	}
	if ( param->ML )
	{
		memory_free(param->ML);
		param->ML = NULL;
	}
	if ( param->MAP )
	{
		memory_free(param->MAP);
		param->MAP = NULL;
	}
	if ( param->prior_rfsize )
	{
		memory_free(param->prior_rfsize);
		param->prior_rfsize = NULL;
	}
	if ( param->prior_rfsize_by_family )
	{
		memory_free_2dim((void**)param->prior_rfsize_by_family, param->pfamily->flist->size, FAMILYSIZEMAX, NULL);
		param->prior_rfsize_by_family = NULL;
	}
	
	int nnodes = param->pcafe ? ((pTree)param->pcafe)->nlist->size : 0;
	viterbi_parameters_clear(&param->viterbi, nnodes);
	if ( !btree_skip && param->pcafe ) 
	{
		if ( param->pcafe->pbdc_array )
		{
			birthdeath_cache_array_free(param->pcafe->pbdc_array);
		}
		cafe_tree_free(param->pcafe);
		//memory_free( param->branchlengths_sorted );
		//param->branchlengths_sorted = NULL;
		memory_free( param->old_branchlength );
		param->old_branchlength = NULL;
		param->pcafe = NULL;
	}
	if ( param->pfamily ) 
	{
		cafe_family_free(param->pfamily);
		param->pfamily = NULL;
	}
	if ( param->parameters ) 
	{
		memory_free( param->parameters );	
		param->parameters = NULL;
	}
	if ( param->lambda ) 
	{
		//memory_free( param->lambda ); param->lambda points to param->parameters
		param->lambda = NULL;
	}
	if ( param->mu ) 
	{
		//memory_free( param->mu );	param->mu points to param->parameters
		param->mu = NULL;
	}
	if ( param->lambda_tree ) 
	{
		phylogeny_free(param->lambda_tree);
		param->lambda_tree = NULL;
	}
	if ( param->mu_tree ) 
	{
		phylogeny_free(param->mu_tree);
		param->mu_tree = NULL;
	}
	
	param->eqbg = 0;
	param->posterior = 0;
    param->num_params = 0;
	param->num_lambdas = 0;
	param->num_mus = 0;
	param->parameterized_k_value = 0;
	param->fixcluster0 = 0;
	param->rootfamily_sizes[0] = 0;
	param->rootfamily_sizes[1] = 1;
	param->family_sizes[0] = 0;
	param->family_sizes[1] = 1;
	param->param_set_func = cafe_shell_set_lambda;
	param->num_threads = 1;
	param->num_random_samples = 1000;
	param->pvalue = 0.01;
	param->bl_augment = 0.5;
}

void cafe_shell_set_sizes()
{
	pCafeTree pcafe = cafe_param->pcafe;
	pcafe->rootfamilysizes[0] = cafe_param->rootfamily_sizes[0];
	pcafe->rootfamilysizes[1] = cafe_param->rootfamily_sizes[1];
	pcafe->familysizes[0] = cafe_param->family_sizes[0];
	pcafe->familysizes[1] = cafe_param->family_sizes[1];
	cafe_tree_set_parameters(pcafe, cafe_param->family_sizes, cafe_param->rootfamily_sizes, 0 );
}

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

void cafe_shell_set_lambda(pCafeParam param, double* parameters)
{
	int i,k;

	if (param->parameters[0] != parameters[0]) memcpy(param->parameters, parameters, param->num_params*sizeof(double));
	// set lambda
	param->lambda = param->parameters;
	// set k_weights
	if (param->parameterized_k_value > 0) {
		double sumofweights = 0;
		for (i = 0; i < (param->parameterized_k_value-1); i++) {
			param->k_weights[i] = param->parameters[param->num_lambdas*(param->parameterized_k_value-param->fixcluster0)+i];
			sumofweights += param->k_weights[i];
		}
		param->k_weights[i] = 1 - sumofweights;
		if( param->p_z_membership == NULL) {
			param->p_z_membership = (double**) memory_new_2dim(param->pfamily->flist->size,param->num_lambdas*param->parameterized_k_value,sizeof(double));
			// assign based on param->k_weights (prior)
			for ( i = 0 ; i < param->pfamily->flist->size ; i++ )
			{
				for (k=0; k<param->parameterized_k_value; k++) {
					param->p_z_membership[i][k] = param->k_weights[k];
				}
			}
		}
	}

	
	param->pcafe->k = param->parameterized_k_value;
	pArrayList nlist = param->pcafe->super.nlist;
	pTree tlambda = param->lambda_tree;
	if ( tlambda == NULL )
	{
		for ( i = 0 ; i < nlist->size ; i++ )
		{
			pCafeNode pcnode = (pCafeNode)nlist->array[i];
			if (param->parameterized_k_value > 0) {
				pcnode->lambda = -1;
				pcnode->mu = -1;
				
				if (pcnode->param_lambdas) { memory_free(pcnode->param_lambdas); pcnode->param_lambdas=NULL;}
				pcnode->param_lambdas = (double*) memory_new(param->parameterized_k_value, sizeof(double));
				if (! param->fixcluster0) {
					memcpy(&pcnode->param_lambdas[0], &parameters[0],(param->parameterized_k_value)*sizeof(double));
				}
				else {
					pcnode->param_lambdas[0] = 0;
					memcpy(&pcnode->param_lambdas[1], &parameters[0],(param->parameterized_k_value-1)*sizeof(double));
				}
				
				//if (pcnode->param_weights) { memory_free(pcnode->param_weights); pcnode->param_weights=NULL;}
				//pcnode->param_weights = (double*) memory_new_with_init(param->k, sizeof(double), param->k_weights);
				
				reset_k_likelihoods(pcnode, param->parameterized_k_value, param->pcafe->size_of_factor);
				
				if (pcnode->k_bd) { arraylist_free(pcnode->k_bd, NULL); }
				pcnode->k_bd = arraylist_new(param->parameterized_k_value);
			}
			else {
				pcnode->lambda = parameters[0];
				pcnode->mu = -1;
			}
		}
	}
	else
	{
		pArrayList lambda_nlist = tlambda->nlist;		
		for ( i = 0 ; i < nlist->size ; i++ )
		{
			pPhylogenyNode pnode = (pPhylogenyNode)lambda_nlist->array[i];
			pCafeNode pcnode = (pCafeNode)nlist->array[i];
			if (param->parameterized_k_value > 0) {
				pcnode->lambda = -1;
				pcnode->mu = -1;
				
				if (pcnode->param_lambdas) { memory_free(pcnode->param_lambdas); pcnode->param_lambdas=NULL;}
				pcnode->param_lambdas = (double*) memory_new(param->parameterized_k_value, sizeof(double));
				if (! param->fixcluster0) {
					memcpy(&pcnode->param_lambdas[0], &parameters[pnode->taxaid*param->parameterized_k_value],(param->parameterized_k_value)*sizeof(double));
				}
				else {
					pcnode->param_lambdas[0] = 0;
					memcpy(&pcnode->param_lambdas[1], &parameters[pnode->taxaid*(param->parameterized_k_value-1)],(param->parameterized_k_value-1)*sizeof(double));
				}
				
				//if (pcnode->param_weights) { memory_free(pcnode->param_weights); pcnode->param_weights=NULL;}
				//pcnode->param_weights = (double*) memory_new_with_init(param->parameterized_k_value, sizeof(double), &param->k_weights[pnode->taxaid*param->parameterized_k_value]);
				
				reset_k_likelihoods(pcnode, param->parameterized_k_value, param->pcafe->size_of_factor);

				if (pcnode->k_bd) { arraylist_free(pcnode->k_bd, NULL); }
				pcnode->k_bd = arraylist_new(param->parameterized_k_value);
			}
			else {
				pcnode->lambda = parameters[pnode->taxaid];
				pcnode->mu = -1;
			}
		}
	}
}



void cafe_shell_set_lambda_mu(pCafeParam param, double* parameters)
{
	int i,k;
	
	if (param->parameters[0] != parameters[0]) {
		memcpy(param->parameters, parameters, param->num_params*sizeof(double));
	}
	// set lambda and mu
	cafe_param->lambda = cafe_param->parameters;
	if (param->parameterized_k_value > 0) {
		cafe_param->mu = &(cafe_param->parameters[param->num_lambdas*(param->parameterized_k_value-param->fixcluster0)]);
	}
	else {
		cafe_param->mu = &(cafe_param->parameters[param->num_lambdas]);
	}
	// set k_weights
	if (param->parameterized_k_value > 0) {
		double sumofweights = 0;
		for (i = 0; i < (param->parameterized_k_value-1); i++) {
			param->k_weights[i] = param->parameters[param->num_lambdas*(param->parameterized_k_value-param->fixcluster0)+(param->num_mus-param->eqbg)*(param->parameterized_k_value-param->fixcluster0)+i];
			sumofweights += param->k_weights[i];
		}
		param->k_weights[i] = 1 - sumofweights;
		if( param->p_z_membership == NULL) {
			param->p_z_membership = (double**) memory_new_2dim(param->pfamily->flist->size,param->num_lambdas*param->parameterized_k_value,sizeof(double));
			// assign based on param->k_weights (prior)
			for ( i = 0 ; i < param->pfamily->flist->size ; i++ )
			{
				for (k=0; k<param->parameterized_k_value; k++) {
					param->p_z_membership[i][k] = param->k_weights[k];
				}
			}
		}
	}
	
	param->pcafe->k = param->parameterized_k_value;
	pArrayList nlist = param->pcafe->super.nlist;
	pTree tlambda = param->lambda_tree;
	if ( tlambda == NULL )
	{
		for ( i = 0 ; i < nlist->size ; i++ )
		{
			pCafeNode pcnode = (pCafeNode)nlist->array[i];
			if (param->parameterized_k_value > 0) {
				pcnode->lambda = -1;
				pcnode->mu = -1;
				
				if (pcnode->param_lambdas) { memory_free(pcnode->param_lambdas); pcnode->param_lambdas=NULL;}
				pcnode->param_lambdas = (double*) memory_new(param->parameterized_k_value, sizeof(double));
				if (! param->fixcluster0) {
					memcpy(&pcnode->param_lambdas[0], &parameters[0],(param->parameterized_k_value)*sizeof(double));
				}
				else {
					pcnode->param_lambdas[0] = 0;
					memcpy(&pcnode->param_lambdas[1], &parameters[0],(param->parameterized_k_value-param->fixcluster0)*sizeof(double));
				}
				
				if (pcnode->param_mus) { memory_free(pcnode->param_mus); pcnode->param_mus = NULL;}
				pcnode->param_mus = (double*) memory_new(param->parameterized_k_value, sizeof(double));
				if (! param->fixcluster0) {
					memcpy(&pcnode->param_mus[0], &parameters[param->num_lambdas*param->parameterized_k_value],(param->parameterized_k_value)*sizeof(double));
				}
				else {
					pcnode->param_mus[0] = 0;
					memcpy(&pcnode->param_mus[1], &parameters[param->num_lambdas*(param->parameterized_k_value-param->fixcluster0)],(param->parameterized_k_value-param->fixcluster0)*sizeof(double));
				}
				
				// (pcnode->param_weights) { memory_free(pcnode->param_weights); pcnode->param_weights=NULL;}
				//pcnode->param_weights = (double*) memory_new_with_init(param->parameterized_k_value, sizeof(double), param->k_weights);
				
				reset_k_likelihoods(pcnode, param->parameterized_k_value, param->pcafe->size_of_factor);

				if (pcnode->k_bd) { arraylist_free(pcnode->k_bd, NULL); }
				pcnode->k_bd = arraylist_new(param->parameterized_k_value);
			}
			else {
				pcnode->lambda = parameters[0];
				pcnode->mu = parameters[param->num_lambdas];
			}
		}
	}
	else
	{
		pArrayList lambda_nlist = tlambda->nlist;		
		for ( i = 0 ; i < nlist->size ; i++ )
		{
			pPhylogenyNode pnode = (pPhylogenyNode)lambda_nlist->array[i];
			pCafeNode pcnode = (pCafeNode)nlist->array[i];
			
			if (param->parameterized_k_value > 0) {
				pcnode->lambda = -1;
				pcnode->mu = -1;
				
				// set lambdas
				if (pcnode->param_lambdas) { memory_free(pcnode->param_lambdas); pcnode->param_lambdas=NULL;}
				pcnode->param_lambdas = (double*) memory_new(param->parameterized_k_value, sizeof(double));
				if (! param->fixcluster0) {
					memcpy(&pcnode->param_lambdas[0], &parameters[pnode->taxaid*param->parameterized_k_value],(param->parameterized_k_value)*sizeof(double));
				}
				else {
					pcnode->param_lambdas[0] = 0;
					memcpy(&pcnode->param_lambdas[1], &parameters[pnode->taxaid*(param->parameterized_k_value-1)],(param->parameterized_k_value-1)*sizeof(double));
				}
				
				// set mus
				if (pcnode->param_mus) { memory_free(pcnode->param_mus); pcnode->param_mus=NULL;}
				pcnode->param_mus = (double*) memory_new(param->parameterized_k_value, sizeof(double));
				if (param->eqbg) {
					if (pnode->taxaid == 0) {
						memcpy(pcnode->param_mus, pcnode->param_lambdas, (param->parameterized_k_value-param->fixcluster0)*sizeof(double));
					}
					else {
						if (! param->fixcluster0) {
							memcpy(&pcnode->param_mus[0], &parameters[(param->num_lambdas)*param->parameterized_k_value+(pnode->taxaid-param->eqbg)*param->parameterized_k_value],(param->parameterized_k_value)*sizeof(double));
						}
						else {
							pcnode->param_mus[0] = 0;
							memcpy(&pcnode->param_mus[1], &parameters[(param->num_lambdas)*(param->parameterized_k_value-1)+(pnode->taxaid-param->eqbg)*(param->parameterized_k_value-1)],(param->parameterized_k_value-1)*sizeof(double));
						}					
					}
				}
				else {
					if (! param->fixcluster0) {
						memcpy(&pcnode->param_mus[0], &parameters[(param->num_lambdas)*param->parameterized_k_value+pnode->taxaid*param->parameterized_k_value],(param->parameterized_k_value)*sizeof(double));
					}
					else {
						pcnode->param_mus[0] = 0;
						memcpy(&pcnode->param_mus[1], &parameters[(param->num_lambdas)*(param->parameterized_k_value-1)+pnode->taxaid*(param->parameterized_k_value-1)],(param->parameterized_k_value-1)*sizeof(double));
					}
				}			
				
				//if (pcnode->param_weights) { memory_free(pcnode->param_weights); pcnode->param_weights=NULL;}
				//pcnode->param_weights = (double*) memory_new_with_init(param->parameterized_k_value, sizeof(double), &param->k_weights[pnode->taxaid*param->parameterized_k_value]);
				
				reset_k_likelihoods(pcnode, param->parameterized_k_value, param->pcafe->size_of_factor);

				if (pcnode->k_bd) { arraylist_free(pcnode->k_bd, NULL); }
				pcnode->k_bd = arraylist_new(param->parameterized_k_value);
			}
			else {
				if (param->eqbg) {
					pcnode->lambda = parameters[pnode->taxaid];
					if (pnode->taxaid == 0) {
						pcnode->mu = pcnode->lambda;
					}
					else {
						pcnode->mu = parameters[(param->num_lambdas)+(pnode->taxaid-param->eqbg)];
					}
				}
				else {
					pcnode->lambda = parameters[pnode->taxaid];
					pcnode->mu = parameters[(param->num_lambdas)+pnode->taxaid];
				}
			}
		}
	}
}


int cafe_shell_parse_familysize(int argc, char* argv[])
{
	int i;
	int max = 0;
	pArrayList nlist = cafe_param->pcafe->super.nlist;
	if ( argc != (nlist->size/2+1) )
	{
		fprintf( stderr, "ERROR: There are %d species\n", nlist->size/2 + 1);
		return -1;
	}
	for ( i = 0; i < argc ; i++ )
	{
		int size = -1;
		sscanf( argv[i],"%d", &size );
		STDERR_IF(size<0, "ERROR: You must input integers greater than or equal to 0\n");
		pCafeNode pnode = (pCafeNode)nlist->array[i*2];
		pnode->familysize = size;
		if ( size > max ) max = size;
	}
	return max;
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

int cafe_shell_parse_branchlength(int argc, char* argv[])
{
	int i;


	pArrayList nlist = cafe_param->pcafe->super.nlist;

	if ( argc != nlist->size )
	{
		fprintf( stderr, "ERROR: There are %d branches including the empty branch of root\n", nlist->size );
		return -1;
	}
	for ( i = 0; i < argc ; i++ )
	{
		int size = -1;
		pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
		sscanf( argv[i],"%d", &size );
		if ( size > 0 )
		{ 
			pnode->branchlength = size;
		}
		else
		{
			fprintf(stderr,"ERROR: the branch length of node %d is not changed\n", i);
		}
	}
	if ( cafe_param->pcafe->pbdc_array ) cafe_tree_set_birthdeath(cafe_param->pcafe);
	return 0;
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
	if ( cafe_param->pcafe->pbdc_array ) cafe_tree_set_birthdeath(cafe_param->pcafe);
	return 0;
}


double* cafe_shell_likelihood(int max)
{
	int rfmax = max+ MAX(50,max/4);
	int fmax = max + MAX(50,max/5);
	pCafeTree pcafe = cafe_param->pcafe;

	if ( pcafe->pbdc_array == NULL || pcafe->pbdc_array->maxFamilysize <  MAX(rfmax, fmax)  )
	{
		cafe_param->rootfamily_sizes[1] = rfmax;
		cafe_param->family_sizes[1] = fmax;
		cafe_tree_set_parameters(pcafe, cafe_param->family_sizes, cafe_param->rootfamily_sizes, 0 );
		if ( pcafe->pbdc_array )
		{
			int remaxFamilysize = MAX(cafe_param->family_sizes[1], cafe_param->rootfamily_sizes[1]);
			birthdeath_cache_resize(cafe_param->pcafe->pbdc_array, remaxFamilysize);
			cafe_tree_set_birthdeath(cafe_param->pcafe);
		}
		else
		{
			cafe_set_birthdeath_cache(cafe_param);
		}
	}
	else
	{
		pcafe->rootfamilysizes[0] = 1;
		pcafe->rootfamilysizes[1] = rfmax;
		pcafe->familysizes[1] = fmax;
		pcafe->rfsize = pcafe->rootfamilysizes[1] - pcafe->rootfamilysizes[0] + 1;
	}
	return cafe_tree_likelihood(pcafe);
}

/**
* \brief Initializes the global \ref cafe_param that holds the data acted upon by cafe. Called at program startup.
* 
*/
void cafe_shell_init(int quiet)
{
	cafe_param = (pCafeParam)memory_new(1,sizeof(CafeParam));
	cafe_param->rootfamily_sizes[0] = 1;
	cafe_param->rootfamily_sizes[1] = 1;
	cafe_param->family_sizes[0] = 0;
	cafe_param->family_sizes[1] = 1;
	cafe_param->param_set_func = cafe_shell_set_lambda;
	cafe_param->flog = stdout;
	cafe_param->num_threads = 1;
	cafe_param->num_random_samples = 1000;
	cafe_param->bl_augment = 0.5;
	cafe_param->pvalue = 0.01;
	cafe_param->quiet = quiet;
	cafe_param->prior_rfsize_by_family = NULL;
}

int cafe_cmd_tree(int argc, char* argv[])
{
    char buf[STRING_BUF_SIZE];
	if ( argc == 1 )
	{
		printf("Newick: ");
		if (fgets(buf,STRING_BUF_SIZE,stdin) == NULL)
			fprintf(stderr, "Failed to read input\n");

	}
	else
	{
		string_pchar_join(buf,NULL,argc-1, &argv[1]);
	}
	if ( cafe_param->pcafe ) 
	{
		if ( cafe_param->pcafe->pbdc_array )
		{
			cafe_free_birthdeath_cache(cafe_param->pcafe);
		}
		cafe_tree_free(cafe_param->pcafe);
		memory_free( cafe_param->old_branchlength );
		cafe_param->old_branchlength = NULL;
	}
	cafe_param->pcafe = cafe_tree_new(buf, cafe_param->family_sizes, 
			                          cafe_param->rootfamily_sizes, 0, 0);
	if (cafe_param->pcafe == NULL) {
		return -1;
	}
	cafe_param->num_branches = cafe_param->pcafe->super.nlist->size - 1;
	cafe_param->old_branchlength = (int*)memory_new(cafe_param->num_branches, sizeof(int));
	pTree ptree = (pTree)cafe_param->pcafe;
	int i=0, j=0;
    // find max_branch_length and sum_branch_length.
    cafe_param->max_branch_length = 0;
    cafe_param->sum_branch_length = 0;
    for( j = 0 ; j < ptree->nlist->size; j++ )
	{
		pPhylogenyNode pnode = (pPhylogenyNode)ptree->nlist->array[j];
		if ( pnode->branchlength > 0 )
		{
            cafe_param->sum_branch_length += pnode->branchlength;
            if (cafe_param->max_branch_length < pnode->branchlength) {
                cafe_param->max_branch_length = pnode->branchlength;
            }
            i++;
		}
	}
	if ( i != ptree->nlist->size - 1 )
	{
		fprintf(stderr, "Expected number of Branch length is %d\n", ptree->nlist->size - 1 );
		cafe_tree_free(cafe_param->pcafe);
		cafe_param->pcafe = NULL;
		return -1;
	}
	if (!cafe_param->quiet)
		cafe_tree_string_print(cafe_param->pcafe);
	if ( cafe_param->pfamily )
	{
		cafe_family_set_species_index(cafe_param->pfamily, cafe_param->pcafe);
	}
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


int cafe_cmd_lambda_mu(int argc, char* argv[])
{
	int i,j;
	STDERR_IF( cafe_param->pfamily == NULL, "ERROR(lambdamu): You must load family data first: command 'load'\n");
	STDERR_IF( cafe_param->pcafe == NULL, "ERROR(lambdamu): You did not specify tree: command 'tree'\n" );
	
	
	pCafeTree pcafe = cafe_param->pcafe;
	pArrayList pargs = cafe_shell_build_argument(argc,argv);
	cafe_param->lambda = NULL;
	cafe_param->mu = NULL;
	if (cafe_param->lambda_tree) {
		phylogeny_free(cafe_param->lambda_tree);
		cafe_param->lambda_tree = NULL;
	}
	if ( cafe_param->mu_tree ) 
	{
		phylogeny_free(cafe_param->mu_tree);
		cafe_param->mu_tree = NULL;
	}
	cafe_param->num_lambdas = -1;
	cafe_param->num_mus = -1;
	cafe_param->parameterized_k_value = 0;
	cafe_param->param_set_func = cafe_shell_set_lambda_mu;
	
	int bdone = 0;
	int bsearch = 0;
	int bprint = 0;
	
	//////
	CafeParam tmpparam;
	pCafeParam tmp_param = &tmpparam;
	memset(tmp_param, 0, sizeof(CafeParam));
    tmp_param->posterior = 1;			
	STDERR_IF( ( cafe_param->pfamily == NULL || cafe_param->pcafe == NULL ), "ERROR(lambda): Please load family (\"load\") and cafe tree (\"tree\") before running \"lambda\" command.");
	
	for ( i = 0 ; i < pargs->size ; i++ )
	{
		pArgument parg = (pArgument)pargs->array[i];
		
		// Search for whole family 
		if ( !strcmp( parg->opt, "-s" ) )
		{
			bsearch = 1;
		}
		else if ( !strcmp( parg->opt, "-checkconv" ) )
		{
			tmp_param->checkconv = 1;
		}		
		else if ( !strcmp( parg->opt, "-t") )
		{
			bdone = __cafe_cmd_lambda_tree(parg);
			if (bdone < 0) {
				return -1;
			}
			pString pstr = phylogeny_string(cafe_param->lambda_tree,NULL);
			cafe_log(cafe_param,"Lambda Tree: %s\n", pstr->buf);
			string_free(pstr);
			tmp_param->lambda_tree = cafe_param->lambda_tree;
			tmp_param->num_lambdas = cafe_param->num_lambdas;
			tmp_param->num_mus = cafe_param->num_mus = cafe_param->num_lambdas;
		}
		else if ( !strcmp( parg->opt, "-l") )
		{
			if ( tmp_param->lambda ) memory_free(tmp_param->lambda);
			tmp_param->lambda = NULL;
			tmp_param->lambda = (double*)memory_new(parg->argc, sizeof(double) );
			for ( j = 0 ; j < parg->argc; j++ )
			{
			 	sscanf( parg->argv[j], "%lf", &tmp_param->lambda[j] );
			}
			tmp_param->num_params += parg->argc;
		}
		else if ( !strcmp( parg->opt, "-m") )
		{
			if ( tmp_param->mu ) memory_free(tmp_param->mu);
			tmp_param->mu = NULL;
			tmp_param->mu = (double*)memory_new(parg->argc, sizeof(double) );
			for ( j = 0 ; j < parg->argc; j++ )
			{
			 	sscanf( parg->argv[j], "%lf", &tmp_param->mu[j] );
			}
			tmp_param->num_params += parg->argc;
		}
		else if ( !strcmp( parg->opt, "-p") )
		{
			if ( tmp_param->k_weights ) memory_free(tmp_param->k_weights);
			tmp_param->k_weights = NULL;
			tmp_param->k_weights = (double*)memory_new(parg->argc, sizeof(double) );
			for ( j = 0 ; j < parg->argc; j++ )
			{
				sscanf( parg->argv[j], "%lf", &tmp_param->k_weights[j] );
			}
			tmp_param->num_params += parg->argc;
		}
		else if ( !strcmp (parg->opt, "-k") ) 
		{
			sscanf( parg->argv[0], "%d", &tmp_param->parameterized_k_value );	
		}
		else if ( !strcmp (parg->opt, "-f") ) 
		{
			tmp_param->fixcluster0 = 1;
		}
		else if ( !strcmp( parg->opt, "-eqbg") ) 
		{
			tmp_param->eqbg = 1;
		}
	}
		
	//////////
	arraylist_free( pargs, free );
	
	if ( bdone ) 
	{
		if ( bdone ) return 0;
	}
	
	// copy parameters collected to cafe_param based on the combination of options.
	{
		cafe_param->posterior = tmp_param->posterior;
		if (cafe_param->posterior) {
			// set rootsize prior based on leaf size
			cafe_set_prior_rfsize_empirical(cafe_param);
		}		
		// search or set
		if (bsearch) {
            // prepare parameters
			if (tmp_param->lambda_tree != NULL) {
				// cafe_param->num_lambdas determined by lambda tree.
				cafe_param->eqbg = tmp_param->eqbg;
				if (tmp_param->parameterized_k_value > 0) {
					cafe_param->parameterized_k_value = tmp_param->parameterized_k_value;
					cafe_param->fixcluster0 = tmp_param->fixcluster0;
					cafe_param->num_params = (tmp_param->num_lambdas*(tmp_param->parameterized_k_value-tmp_param->fixcluster0))+
					((tmp_param->num_mus-tmp_param->eqbg)*(tmp_param->parameterized_k_value-tmp_param->fixcluster0))+
					(tmp_param->parameterized_k_value-1);
					
					if( cafe_param->parameters ) memory_free(cafe_param->parameters);
					cafe_param->parameters = NULL;
					cafe_param->parameters = (double*)memory_new(cafe_param->num_params, sizeof(double));
					if (cafe_param->k_weights) { memory_free(cafe_param->k_weights);}
					cafe_param->k_weights = NULL;
					cafe_param->k_weights = (double*) memory_new(cafe_param->parameterized_k_value, sizeof(double));
				}
				else {	// search whole dataset branch specific
					cafe_param->num_params = tmp_param->num_lambdas+(tmp_param->num_mus-tmp_param->eqbg);
					
					if( cafe_param->parameters ) memory_free(cafe_param->parameters);
					cafe_param->parameters = NULL;
					cafe_param->parameters = (double*)memory_new(cafe_param->num_params, sizeof(double));
				}
			}
			else {
				cafe_param->num_lambdas = tmp_param->num_lambdas = 1;	
				cafe_param->num_mus = tmp_param->num_mus = 1;	
				if (tmp_param->eqbg) {
					fprintf( stderr, "ERROR(lambdamu): Cannot use option eqbg without specifying a lambda tree. \n");
					return -1;											
				}
				if (tmp_param->parameterized_k_value > 0) {
					cafe_param->parameterized_k_value = tmp_param->parameterized_k_value;
					cafe_param->fixcluster0 = tmp_param->fixcluster0;
					cafe_param->num_params = (tmp_param->num_lambdas*(tmp_param->parameterized_k_value-tmp_param->fixcluster0))+
					(tmp_param->num_mus*(tmp_param->parameterized_k_value-tmp_param->fixcluster0))+
					(tmp_param->parameterized_k_value-1);
					
					if( cafe_param->parameters ) memory_free(cafe_param->parameters);
					cafe_param->parameters = NULL;
					cafe_param->parameters = (double*)memory_new(cafe_param->num_params, sizeof(double));
					if (cafe_param->k_weights) { memory_free(cafe_param->k_weights);}
					cafe_param->k_weights = NULL;
					cafe_param->k_weights = (double*) memory_new(cafe_param->parameterized_k_value, sizeof(double));
				}
				else {	// search whole dataset whole tree
					cafe_param->num_params = tmp_param->num_lambdas+tmp_param->num_mus;
					
					if( cafe_param->parameters ) memory_free(cafe_param->parameters);
					cafe_param->parameters = NULL;
					cafe_param->parameters = (double*)memory_new(cafe_param->num_params, sizeof(double));
				}
			}
			// search
			if (tmp_param->checkconv) { cafe_param->checkconv = 1; }
			cafe_best_lambda_mu_by_fminsearch(cafe_param, cafe_param->num_lambdas, cafe_param->num_mus, cafe_param->parameterized_k_value);
		}
		else {
			if (tmp_param->lambda_tree != NULL) {
				// cafe_param->num_lambdas determined by lambda tree.
				cafe_param->eqbg = tmp_param->eqbg;
				if (tmp_param->parameterized_k_value > 0) {	// search clustered branch specific
					cafe_param->parameterized_k_value = tmp_param->parameterized_k_value;
					cafe_param->fixcluster0 = tmp_param->fixcluster0;
					cafe_param->num_params = (tmp_param->num_lambdas*(tmp_param->parameterized_k_value-tmp_param->fixcluster0))+
					((tmp_param->num_mus-tmp_param->eqbg)*(tmp_param->parameterized_k_value-tmp_param->fixcluster0))+
					(tmp_param->parameterized_k_value-1);
					
					// check if the numbers of lambdas and proportions put in matches the number of parameters
					if (cafe_param->num_params != tmp_param->num_params) {
						fprintf( stderr, "ERROR(lambdamu): Number of parameters not correct. \n");
						fprintf( stderr, "the number of -l lambdas -m mus and -p proportions are %d they need to be %d\n", tmp_param->num_params, cafe_param->num_params );
						pString pstr = phylogeny_string(tmp_param->lambda_tree,NULL);
						fprintf( stderr, "based on the tree %s and -k clusters %d.\n", pstr->buf, cafe_param->parameterized_k_value );
						string_free(pstr);
						return -1;						
					}
					
					// copy user input into parameters
					if( cafe_param->parameters ) memory_free(cafe_param->parameters);
					cafe_param->parameters = NULL;
					cafe_param->parameters = (double*)memory_new(cafe_param->num_params, sizeof(double));
					memcpy(cafe_param->parameters,tmp_param->lambda, sizeof(double)*tmp_param->num_lambdas*(tmp_param->parameterized_k_value-tmp_param->fixcluster0));
					memcpy(&cafe_param->parameters[cafe_param->num_lambdas*(cafe_param->parameterized_k_value-tmp_param->fixcluster0)],tmp_param->mu, sizeof(double)*((tmp_param->num_mus-tmp_param->eqbg)*(tmp_param->parameterized_k_value-tmp_param->fixcluster0)));
					memcpy(&cafe_param->parameters[(cafe_param->num_lambdas*(cafe_param->parameterized_k_value-tmp_param->fixcluster0))+((tmp_param->num_mus-tmp_param->eqbg)*(tmp_param->parameterized_k_value-tmp_param->fixcluster0))], tmp_param->k_weights, sizeof(double)*(tmp_param->parameterized_k_value-1));
					// prepare space for k_weights
					if ( cafe_param->k_weights ) memory_free(cafe_param->k_weights);
					cafe_param->k_weights = NULL;
					cafe_param->k_weights = (double*)memory_new(cafe_param->parameterized_k_value-1, sizeof(double) );										
				}
				else {	// search whole dataset branch specific
					cafe_param->num_params = tmp_param->num_lambdas+(tmp_param->num_mus-tmp_param->eqbg);
					
					// check if the numbers of lambdas and proportions put in matches the number of parameters
					if (cafe_param->num_params != tmp_param->num_params) {
						fprintf( stderr, "ERROR(lambdamu): Number of parameters not correct. \n");
						fprintf( stderr, "the number of -l lambdas -m mus are %d they need to be %d\n", tmp_param->num_params, cafe_param->num_params );
						pString pstr = phylogeny_string(tmp_param->lambda_tree,NULL);
						fprintf( stderr, "based on the tree %s \n", pstr->buf );
						string_free(pstr);
						return -1;						
					}
					
					// copy user input into parameters
					if( cafe_param->parameters ) memory_free(cafe_param->parameters);
					cafe_param->parameters = NULL;
					cafe_param->parameters = (double*)memory_new(cafe_param->num_params, sizeof(double));
					memcpy(cafe_param->parameters,tmp_param->lambda, sizeof(double)*tmp_param->num_lambdas);
					memcpy(&cafe_param->parameters[cafe_param->num_lambdas],tmp_param->mu, sizeof(double)*(tmp_param->num_mus-tmp_param->eqbg));
					
				}
			}
			else {
				cafe_param->num_lambdas = tmp_param->num_lambdas = 1;	
				cafe_param->num_mus = tmp_param->num_mus = 1;	
				if (tmp_param->eqbg) {
					fprintf( stderr, "ERROR(lambdamu): Cannot use option eqbg without specifying a lambda tree. \n");
					return -1;											
				}
				if (tmp_param->parameterized_k_value > 0) {				// search clustered whole tree
					cafe_param->parameterized_k_value = tmp_param->parameterized_k_value;
					cafe_param->fixcluster0 = tmp_param->fixcluster0;
					cafe_param->num_params = (tmp_param->num_lambdas*(tmp_param->parameterized_k_value-tmp_param->fixcluster0))+
					(tmp_param->num_mus*(tmp_param->parameterized_k_value-tmp_param->fixcluster0))+
					(tmp_param->parameterized_k_value-1);
					
					// check if the numbers of lambdas and proportions put in matches the number of parameters
					if (cafe_param->num_params != tmp_param->num_params) {
						fprintf( stderr, "ERROR(lambdamu): Number of parameters not correct. \n");
						fprintf( stderr, "the number of -l lambdas -m mus and -p proportions are %d they need to be %d\n", tmp_param->num_params, cafe_param->num_params );
						fprintf( stderr, "based on the -k clusters %d.\n", cafe_param->parameterized_k_value );
						return -1;						
					}
					
					// copy user input into parameters
					if( cafe_param->parameters ) memory_free(cafe_param->parameters);
					cafe_param->parameters = NULL;
					cafe_param->parameters = (double*)memory_new(cafe_param->num_params, sizeof(double));
					memcpy(cafe_param->parameters,tmp_param->lambda, sizeof(double)*tmp_param->num_lambdas*(tmp_param->parameterized_k_value-tmp_param->fixcluster0));
					memcpy(&cafe_param->parameters[cafe_param->num_lambdas*(cafe_param->parameterized_k_value-tmp_param->fixcluster0)], tmp_param->mu, sizeof(double)*tmp_param->num_mus*(cafe_param->parameterized_k_value-tmp_param->fixcluster0));
					memcpy(&cafe_param->parameters[cafe_param->num_lambdas*(cafe_param->parameterized_k_value-tmp_param->fixcluster0)+tmp_param->num_mus*(cafe_param->parameterized_k_value-tmp_param->fixcluster0)], tmp_param->k_weights, sizeof(double)*(tmp_param->parameterized_k_value-1));
					// prepare space for k_weights
					if ( cafe_param->k_weights ) memory_free(cafe_param->k_weights);
					cafe_param->k_weights = NULL;
					cafe_param->k_weights = (double*)memory_new(cafe_param->parameterized_k_value-1, sizeof(double) );										
					
				}
				else {	// search whole dataset whole tree
					cafe_param->num_params = tmp_param->num_lambdas+tmp_param->num_mus;
					
					// check if the numbers of lambdas and proportions put in matches the number of parameters
					if (cafe_param->num_params != tmp_param->num_params) {
						fprintf( stderr, "ERROR(lambdamu): Number of parameters not correct. \n");
						fprintf( stderr, "the number of -l lambdas -m mus are %d they need to be %d\n", tmp_param->num_params, cafe_param->num_params );
						return -1;						
					}
					
					// copy user input into parameters
					if( cafe_param->parameters ) memory_free(cafe_param->parameters);
					cafe_param->parameters = NULL;
					cafe_param->parameters = (double*)memory_new(cafe_param->num_params, sizeof(double));
					memcpy(cafe_param->parameters,tmp_param->lambda, sizeof(double)*tmp_param->num_lambdas);
					memcpy(&cafe_param->parameters[tmp_param->num_lambdas],tmp_param->mu, sizeof(double)*tmp_param->num_mus);
				}
			}
			cafe_param->param_set_func(cafe_param, cafe_param->parameters);
		}
	}
	
	
	//////////
	
	
	if ( bprint )
	{
		pString pstr = cafe_tree_string_with_lambda(pcafe);
		printf("%s\n", pstr->buf );
		string_free(pstr);
	}
	  if ( cafe_param->pfamily )
	 {
	 cafe_set_birthdeath_cache_thread(cafe_param->pcafe, cafe_param->parameterized_k_value, cafe_param->family_sizes, cafe_param->rootfamily_sizes);
	 }
	   
	cafe_log(cafe_param,"DONE: Lamda,Mu Search or setting, for command:\n");
	char buf[STRING_STEP_SIZE];
	buf[0] = '\0';
	string_pchar_join(buf," ",argc, argv);
	cafe_log(cafe_param,"%s\n", buf);
	
	if (bsearch && (cafe_param->parameterized_k_value > 0)) {
		// print the cluster memberships
		cafe_family_print_cluster_membership(cafe_param);
	}
	return 0;
}






void __cafe_cmd_viterbi_print(int max)
{
	pCafeTree pcafe = cafe_param->pcafe;
	double* lh = cafe_shell_likelihood(max);
	double mlh =  __max(lh,pcafe->rfsize);
	cafe_tree_viterbi(pcafe);
	pString pstr = cafe_tree_string(pcafe);
	printf("%g\t%s\n", mlh , pstr->buf );
	string_free(pstr);
}

void __cafe_cmd_viterbi_family_print(int idx)
{
	pCafeTree pcafe = cafe_param->pcafe;
	cafe_family_set_size_with_family_forced(cafe_param->pfamily,idx,pcafe);
	cafe_tree_likelihood(pcafe);
	int ridx =  __maxidx(((pCafeNode)pcafe->super.root)->likelihoods,pcafe->rfsize) + pcafe->rootfamilysizes[0];
	double mlh =  __max( ((pCafeNode)pcafe->super.root)->likelihoods,pcafe->rfsize);
	//cafe_tree_likelihood(pcafe);
	cafe_tree_viterbi(pcafe);
	pString pstr = cafe_tree_string(pcafe);
	printf("%g(%d)\t%s\n", mlh , ridx,  pstr->buf );
	string_free(pstr);
}

int cafe_cmd_viterbi(int argc, char* argv[])
{
	int i;
	pCafeTree pcafe = cafe_param->pcafe;
	STDERR_IF( cafe_param->pcafe == NULL, "ERROR(viterbi): You did not specify tree: command 'tree'\n" );
	STDERR_IF( cafe_param->lambda == NULL, "ERROR(viterbi): You did not estimate the parameters: use command 'lambda'\n" );
	pArrayList pargs = cafe_shell_build_argument(argc, argv);
	pArgument parg;

	if ( (parg = cafe_shell_get_argument( "-all", pargs) )  )
	{
		cafe_shell_set_sizes();
		STDERR_IF( cafe_param->pfamily == NULL, "ERROR(viterbi): You must load family data first: command 'load'\n");
		FILE *fp = stdout;
		if  ( parg->argc )
		{
			fp = (FILE*)fopen(parg->argv[0], "w");
			if ( fp == NULL )
			{
				fprintf(stderr,"ERROR(viterbi): Cannot open %s in write mode.\n", parg->argv[0] );
				return -1;
			}
		}
		pCafeFamilyItem pitem;
/*
		int j;
	   	int	free_v;
		int nnodes = (pcafe->super.nlist->size-1)/2;
		int nrows = cafe_param->pfamily->flist->size ;
		if ( cafe_param->viterbiNodeFamilysizes == NULL )
		{
			free_v = 1;
			cafe_param->viterbiNodeFamilysizes 
				= (int**)memory_new_2dim(nnodes,nrows,sizeof(int));
		}
		else
		{
			free_v = 0;
		}
*/
		double score = 0;
		for ( i = 0 ; i < cafe_param->pfamily->flist->size ; i++ )
		{
			pitem = (pCafeFamilyItem)cafe_param->pfamily->flist->array[i];
			pitem->maxlh = -1;
			cafe_family_set_size_with_family(cafe_param->pfamily,i, pcafe);
			cafe_tree_likelihood(pcafe);
			int ridx =  __maxidx(((pCafeNode)pcafe->super.root)->likelihoods,pcafe->rfsize) + pcafe->rootfamilysizes[0];
			double mlh =  __max( ((pCafeNode)pcafe->super.root)->likelihoods,pcafe->rfsize);
			score += log(mlh);
			cafe_tree_viterbi(pcafe);
/*
			for ( j = 0 ; j < nnodes ; j++ )
			{
				pCafeNode pcnode = (pCafeNode)pcafe->super.nlist->array[2*j+1];
				cafe_param->viterbiNodeFamilysizes[j][i] = pcnode->familysize;
			}
*/
			pString pstr = cafe_tree_string(pcafe);
			fprintf(fp,"%s\t%g\t%s\t%d\n", pitem->id,  mlh, pstr->buf, ridx );
			string_free(pstr);
		}
		cafe_log( cafe_param, "Score: %f\n", score );
		if ( fp != stdout )
		{
			fprintf(fp,"Score: %f\n", score );
		}
/*
		if ( free_v )
		{
			memory_free_2dim((void**)cafe_param->viterbiNodeFamilysizes,nnodes, 0, NULL );
			cafe_param->viterbiNodeFamilysizes = NULL;
		}
*/
		if ( fp != stdout ) fclose( fp );
	}
	else if ( (parg = cafe_shell_get_argument( "-id", pargs) )  )
	{
		int idx = cafe_family_get_index(cafe_param->pfamily, parg->argv[0] );
		if ( idx == -1 )
		{
			fprintf(stderr, "ERROR(viterbi): %s not found\n", parg->argv[0] );
			return -1;
		}
		__cafe_cmd_viterbi_family_print(idx);
	}
	else if ( (parg = cafe_shell_get_argument( "-idx", pargs) )  )
	{
		int idx = -1; sscanf( parg->argv[0], "%d", &idx );
		if ( idx == -1 )
		{
			fprintf(stderr,"ERROR(viterbi): It is not integer\n");
			return -1;
		}
		if ( idx > cafe_param->pfamily->flist->size )
		{
			fprintf(stderr, "ERROR(viterbi): Out of range[0~%d]: %s.\n", cafe_param->pfamily->flist->size, parg->argv[0] ); 
			return -1;
		}
		__cafe_cmd_viterbi_family_print(idx);
	}
	else if ( argc == 1 )
	{
		__cafe_cmd_viterbi_print( cafe_shell_set_familysize() );
	}
	else
	{
		__cafe_cmd_viterbi_print( cafe_shell_parse_familysize( argc-1, &argv[1]) );
	}
	arraylist_free(pargs, free);
	return 0;
}

int cafe_cmd_reconstruction_accuracy(int argc, char* argv[])
{
	int i, j;
	pCafeTree pcafe = cafe_param->pcafe;
	pArrayList pargs = cafe_shell_build_argument(argc, argv);
	pArgument parg;
	double SSE = 0;
	double MSE = 0;
	int nodecnt = 0;
	int errnodecnt = 0;
	if ((parg = cafe_shell_get_argument("-i", pargs)))
	{
		// read in truth data
		pString truthfile = string_join(" ",parg->argc, parg->argv );
		pCafeFamily truthfamily = cafe_family_new( truthfile->buf, 1 );
		if ( truthfamily == NULL ) {
			fprintf(stderr, "failed to read in true values %s\n", truthfile->buf);
			return -1;
		}
		pCafeTree truthtree = cafe_tree_copy(pcafe);
		// set parameters
		if ( truthtree )
		{
			cafe_family_set_species_index(truthfamily, truthtree);
            // compare inferred vs. truth
            for(i=0; i< cafe_param->pfamily->flist->size; i++) 
            {
                cafe_family_set_truesize_with_family(truthfamily,i, truthtree);
                cafe_family_set_size(cafe_param->pfamily,i, pcafe);
                if (cafe_param->posterior) {
                    cafe_tree_viterbi_posterior(pcafe, cafe_param);
                }
                else {
                    cafe_tree_viterbi(pcafe);
                }
                // inner node SSE
                for (j=1; j<pcafe->super.nlist->size; j=j+2) {
                    int error = ((pCafeNode)truthtree->super.nlist->array[j])->familysize-((pCafeNode)pcafe->super.nlist->array[j])->familysize;
                    SSE += pow(error, 2);
                    if (error != 0) { errnodecnt++; }  
                    nodecnt++;
                }
            }
            MSE = SSE/nodecnt;
        }
		
	}
	arraylist_free(pargs, free);
	cafe_log( cafe_param, "ancestral reconstruction SSE %f MSE %f totalnodes %d errornodes %d\n", SSE, MSE, nodecnt, errnodecnt );
	return 0;
}




double _cafe_cross_validate_by_family(char* queryfile, char* truthfile, char* errortype) 
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
	cafe_set_birthdeath_cache_thread(cafe_param->pcafe, cafe_param->parameterized_k_value, cafe_param->family_sizes, cafe_param->rootfamily_sizes);
	
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




double _cafe_cross_validate_by_species(char* validatefile, char* errortype) 
{
	int i, j;
	cafe_family_read_validate_species( cafe_param, validatefile );
	if ( cafe_param->cv_test_count_list == NULL ) return -1;
	// now compare reconstructed count to true count	
	pCafeTree pcafe = cafe_param->pcafe;
	cafe_set_birthdeath_cache_thread(cafe_param->pcafe, cafe_param->parameterized_k_value, cafe_param->family_sizes, cafe_param->rootfamily_sizes);
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







int cafe_cmd_crossvalidation_by_family(int argc, char* argv[])
{
	int i;
	STDERR_IF( cafe_param->pfamily == NULL, "ERROR(cvfamily): You did not load family: command 'load'\n" );
	STDERR_IF( cafe_param->pcafe == NULL, "ERROR(cvfamily): You did not specify tree: command 'tree'\n" );
	STDERR_IF( cafe_param->lambda == NULL, "ERROR(cvfamily): You did not set the parameters: command 'lambda' or 'lambdamu'\n" );
	
	double MSE_allfolds = 0;
	pCafeFamily pcafe_original = cafe_param->pfamily;
	
	if ( argc < 2 )
	{
		fprintf(stderr, "Usage(cvfamily): %s -fold <num>\n", argv[0] );
		return -1;
	}
	pArrayList pargs = cafe_shell_build_argument(argc, argv);
	pArgument parg;
	if ((parg = cafe_shell_get_argument("-fold", pargs)))
	{
		sscanf( parg->argv[0], "%d", &cafe_param->cv_fold );
	}
	// set up the training-validation set
	cafe_family_split_cvfiles_byfamily(cafe_param);
	
	//
	for (i=0; i<cafe_param->cv_fold; i++) 
	{
		char trainfile[STRING_BUF_SIZE];
		char queryfile[STRING_BUF_SIZE];
		char validatefile[STRING_BUF_SIZE];
		sprintf(trainfile, "%s.%d.train", cafe_param->str_fdata->buf, i+1);
		sprintf(queryfile, "%s.%d.query", cafe_param->str_fdata->buf, i+1);
		sprintf(validatefile, "%s.%d.valid", cafe_param->str_fdata->buf, i+1);
		
		// read in training data
		pCafeFamily tmpfamily = cafe_family_new( trainfile, 1 );
		if ( tmpfamily == NULL ) {
			fprintf(stderr, "failed to read in training data %s\n", trainfile);
			return -1;
		}
		cafe_param->pfamily = tmpfamily;
		
		cafe_param->rootfamily_sizes[0] = 1;
		cafe_param->rootfamily_sizes[1] =  MAX(30,rint(cafe_param->pfamily->max_size * 1.25));
		cafe_param->family_sizes[1] = cafe_param->pfamily->max_size + MAX(50,cafe_param->pfamily->max_size/5);
		if ( cafe_param->pcafe )
		{
			cafe_tree_set_parameters(cafe_param->pcafe, cafe_param->family_sizes, cafe_param->rootfamily_sizes, 0 );
			cafe_family_set_species_index(cafe_param->pfamily, cafe_param->pcafe);
		}
		// re-train 
		if (cafe_param->num_mus > 0) {
			cafe_best_lambda_mu_by_fminsearch(cafe_param, cafe_param->num_lambdas, cafe_param->num_mus, cafe_param->parameterized_k_value);
		}
		else {
			cafe_best_lambda_by_fminsearch(cafe_param, cafe_param->num_lambdas, cafe_param->parameterized_k_value);
		}
		
		//cross-validate
		double MSE = _cafe_cross_validate_by_family(queryfile, validatefile, "MSE");
		MSE_allfolds += MSE;
		cafe_log( cafe_param, "MSE fold %d %f\n", i+1, MSE );
		
		cafe_family_free(tmpfamily);
	}
	MSE_allfolds = MSE_allfolds/cafe_param->cv_fold;
	cafe_log( cafe_param, "MSE all folds %f\n", MSE_allfolds);
	
	//re-load the original family file
	cafe_param->pfamily = pcafe_original;
	cafe_param->rootfamily_sizes[0] = 1;
	cafe_param->rootfamily_sizes[1] =  MAX(30,rint(cafe_param->pfamily->max_size * 1.25));
	cafe_param->family_sizes[1] = cafe_param->pfamily->max_size + MAX(50,cafe_param->pfamily->max_size/5);
	if ( cafe_param->pcafe )
	{
		cafe_tree_set_parameters(cafe_param->pcafe, cafe_param->family_sizes, cafe_param->rootfamily_sizes, 0 );
		cafe_family_set_species_index(cafe_param->pfamily, cafe_param->pcafe);
	}
	// re-train 
	if (cafe_param->num_mus > 0) {
		cafe_best_lambda_mu_by_fminsearch(cafe_param, cafe_param->num_lambdas, cafe_param->num_mus, cafe_param->parameterized_k_value);
	}
	else {
		cafe_best_lambda_by_fminsearch(cafe_param, cafe_param->num_lambdas, cafe_param->parameterized_k_value);
	}
	
	// remove training-validation set
	cafe_family_clean_cvfiles_byfamily(cafe_param);	
	return 0;
}


int cafe_cmd_crossvalidation_by_species(int argc, char* argv[])
{
	int i;
	STDERR_IF( cafe_param->pfamily == NULL, "ERROR(cvspecies): You did not load family: command 'load'\n" );
	STDERR_IF( cafe_param->pcafe == NULL, "ERROR(cvspecies): You did not specify tree: command 'tree'\n" );
	STDERR_IF( cafe_param->lambda == NULL, "ERROR(cvspecies): You did not set the parameters: command 'lambda' or 'lambdamu'\n" );
	
	double MSE_allspecies = 0;
	pCafeFamily pcafe_original = cafe_param->pfamily;
	int num_species_original = cafe_param->pfamily->num_species;
	char** species_names_original = cafe_param->pfamily->species;
	
	if ( argc < 2 )
	{
		// set up the training-validation set
		cafe_family_split_cvfiles_byspecies(cafe_param);
		
		for (i=0; i<num_species_original; i++) {
			char trainfile[STRING_BUF_SIZE];
			char validatefile[STRING_BUF_SIZE];
			sprintf(trainfile, "%s.%s.train", cafe_param->str_fdata->buf, species_names_original[i]);
			sprintf(validatefile, "%s.%s.valid", cafe_param->str_fdata->buf, species_names_original[i]);
			
			// read in training data
			pCafeFamily tmpfamily = cafe_family_new( trainfile, 1 );
			if ( tmpfamily == NULL ) {
				fprintf(stderr, "failed to read in training data %s\n", trainfile);
				fprintf(stderr, "did you load the family data with the cross-validation option (load -i <familyfile> -cv)?\n");
				return -1;
			}
			cafe_param->pfamily = tmpfamily;
			
			cafe_param->rootfamily_sizes[0] = 1;
			cafe_param->rootfamily_sizes[1] =  MAX(30,rint(cafe_param->pfamily->max_size * 1.25));
			cafe_param->family_sizes[1] = cafe_param->pfamily->max_size + MAX(50,cafe_param->pfamily->max_size/5);
			if ( cafe_param->pcafe )
			{
				cafe_tree_set_parameters(cafe_param->pcafe, cafe_param->family_sizes, cafe_param->rootfamily_sizes, 0 );
				cafe_family_set_species_index(cafe_param->pfamily, cafe_param->pcafe);
			}
			// re-train 
			if (cafe_param->num_mus > 0) {
				cafe_best_lambda_mu_by_fminsearch(cafe_param, cafe_param->num_lambdas, cafe_param->num_mus, cafe_param->parameterized_k_value);
			}
			else {
				cafe_best_lambda_by_fminsearch(cafe_param, cafe_param->num_lambdas, cafe_param->parameterized_k_value);
			}
			
			//cross-validate
			double MSE = _cafe_cross_validate_by_species(validatefile, "MSE");
			MSE_allspecies += MSE;
			cafe_log( cafe_param, "MSE %s %f\n", cafe_param->cv_species_name, MSE );
			
			cafe_family_free(tmpfamily);
		}
		MSE_allspecies = MSE_allspecies/num_species_original;
		cafe_log( cafe_param, "MSE all species %f\n", MSE_allspecies );
		
		//re-load the original family file
		cafe_param->pfamily = pcafe_original;
		cafe_param->rootfamily_sizes[0] = 1;
		cafe_param->rootfamily_sizes[1] =  MAX(30,rint(cafe_param->pfamily->max_size * 1.25));
		cafe_param->family_sizes[1] = cafe_param->pfamily->max_size + MAX(50,cafe_param->pfamily->max_size/5);
		if ( cafe_param->pcafe )
		{
			cafe_tree_set_parameters(cafe_param->pcafe, cafe_param->family_sizes, cafe_param->rootfamily_sizes, 0 );
			cafe_family_set_species_index(cafe_param->pfamily, cafe_param->pcafe);
		}
		// re-train 
		if (cafe_param->num_mus > 0) {
			cafe_best_lambda_mu_by_fminsearch(cafe_param, cafe_param->num_lambdas, cafe_param->num_mus, cafe_param->parameterized_k_value);
		}
		else {
			cafe_best_lambda_by_fminsearch(cafe_param, cafe_param->num_lambdas, cafe_param->parameterized_k_value);
		}
		
		// remove training-validation set
		cafe_family_clean_cvfiles_byspecies(cafe_param);
	}
	else 
	{
		pArrayList pargs = cafe_shell_build_argument(argc, argv);
		pArgument parg;
		pString str_fdata;
		if ((parg = cafe_shell_get_argument("-i", pargs)))
		{
			str_fdata = string_join(" ",parg->argc, parg->argv );
			_cafe_cross_validate_by_species(str_fdata->buf, "MSE"); 
		}
	}		
	return 0;
}



int cafe_cmd_load(int argc, char* argv[])
{
	if ( argc < 2 )
	{
		fprintf(stderr, "Usage(load): %s <famliy file>\n", argv[0] );
		return -1;
	}
	cafe_shell_clear_param(cafe_param, 1);
	int bfilter = 0;

	pArrayList pargs = cafe_shell_build_argument(argc, argv);
	pArgument parg;

	if ((parg = cafe_shell_get_argument("-t", pargs)))
	{
		sscanf( parg->argv[0], "%d", &cafe_param->num_threads );
	}
	if ((parg = cafe_shell_get_argument("-r", pargs)))
	{
		sscanf( parg->argv[0], "%d", &cafe_param->num_random_samples );
	}
	if ((parg = cafe_shell_get_argument("-p", pargs)))
	{
		sscanf( parg->argv[0], "%lf", &cafe_param->pvalue );
	}
	if ((parg = cafe_shell_get_argument("-l", pargs)))
	{
		__cafe_cmd_log( parg->argc, parg->argv );
	}
	if ((parg = cafe_shell_get_argument("-filter", pargs)))
	{
		if ( cafe_param->pcafe == NULL )
		{
			fprintf(stderr,"Error(load): You did not specify tree. Skip filtering\n");
		}
		bfilter = 1;
	}
	if ((parg = cafe_shell_get_argument("-i", pargs)))
	{
		cafe_param->str_fdata = string_join(" ",parg->argc, parg->argv );
		cafe_param->pfamily = cafe_family_new( cafe_param->str_fdata->buf, 1 );
		if ( cafe_param->pfamily == NULL ) return -1;
	}
	arraylist_free(pargs, free);

	if ( cafe_param->pfamily == NULL )
	{
		fprintf(stderr,"ERROR(load): You must use -i option for input file\n");
		cafe_shell_clear_param(cafe_param, 1);
		return -1;
	}
	
	// this is very important!!!!!!!!!!!
	cafe_param->rootfamily_sizes[0] = 1;
	//cafe_param->rootfamily_sizes[0] = 0;
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	cafe_param->rootfamily_sizes[1] =  MAX(30,rint(cafe_param->pfamily->max_size * 1.25));
	cafe_param->family_sizes[1] = cafe_param->pfamily->max_size + MAX(50,cafe_param->pfamily->max_size/5);
	if ( cafe_param->pcafe )
	{
		cafe_tree_set_parameters(cafe_param->pcafe, cafe_param->family_sizes, cafe_param->rootfamily_sizes, 0 );
		cafe_family_set_species_index(cafe_param->pfamily, cafe_param->pcafe);
	}
	cafe_param->ML = (double*)memory_new( cafe_param->pfamily->flist->size, sizeof(double));
	cafe_param->MAP = (double*)memory_new( cafe_param->pfamily->flist->size, sizeof(double));
	
	if ( cafe_param->pcafe && bfilter )
	{
		cafe_family_filter(cafe_param);
	}
	if ( cafe_param->lambda && cafe_param->pcafe )
	{
		cafe_set_birthdeath_cache_thread(cafe_param->pcafe, cafe_param->parameterized_k_value, cafe_param->family_sizes, cafe_param->rootfamily_sizes);
	}
	cafe_cmd_print_param(0,NULL);
	return 0;
}

int cafe_cmd_save(int argc, char* argv[])
{
	if ( argc != 2 )
	{
		fprintf(stderr,"Usage(save): save filename\n");
		return -1;
	}
	FILE* fp;
	if ( ( fp = fopen( argv[1], "w" ) ) == NULL )
	{
		fprintf(stderr,"ERROR(save): Cannot open %s in write mode.\n", argv[1] );
		return -1;
	}
	pCafeFamily pcf = cafe_param->pfamily;

	int i, n;
	char buf[STRING_STEP_SIZE];
	string_pchar_join(buf,"\t", pcf->num_species, pcf->species );
	fprintf(fp,"Desc\tFamily ID\t%s\n", buf );	

	for ( i = 0 ; i < pcf->flist->size ; i++ )
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[i];
		fprintf(fp,"%s\t%s\t%d", pitem->desc,pitem->id, pitem->count[0]);	
		for ( n =  1 ; n < pcf->num_species ; n++ )
		{
			fprintf(fp,"\t%d", pitem->count[n]);	
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
	return 0;
}

int cafe_cmd_print_param(int argc, char* argv[])
{
	cafe_log( cafe_param, "-----------------------------------------------------------\n");
	cafe_log( cafe_param, "Family information: %s\n", cafe_param->str_fdata->buf );
	cafe_log( cafe_param, "Log: %s\n", cafe_param->flog == stdout ? "stdout" : cafe_param->str_log->buf );
	if( cafe_param->pcafe ) 
	{
		pString pstr = phylogeny_string( (pTree)cafe_param->pcafe, NULL );
		cafe_log( cafe_param, "Tree: %s\n", pstr->buf ); 
		string_free(pstr);
	}
	cafe_log( cafe_param, "The number of families is %d\n", cafe_param->pfamily->flist->size );
	cafe_log( cafe_param, "Root Family size : %d ~ %d\n", cafe_param->rootfamily_sizes[0], cafe_param->rootfamily_sizes[1]);
	cafe_log( cafe_param, "Family size : %d ~ %d\n", cafe_param->family_sizes[0], cafe_param->family_sizes[1]);
	cafe_log( cafe_param, "P-value: %f\n", cafe_param->pvalue);
	cafe_log( cafe_param, "Num of Threads: %d\n", cafe_param->num_threads );
	cafe_log( cafe_param, "Num of Random: %d\n", cafe_param->num_random_samples );
	if ( cafe_param->lambda )
	{
		pString pstr = cafe_tree_string_with_lambda(cafe_param->pcafe);
		cafe_log( cafe_param, "Lambda: %s\n", pstr->buf );
		string_free(pstr);
	}
	return 0;
}

int cafe_cmd_version(int argc, char* argv[])
{
	printf("Version: %g, at %s\n", 3.1, __DATE__ );
	return 0;
}

void __cafe_cmd_pvalue_print(int max)
{
	pCafeTree pcafe = cafe_param->pcafe;
	double* lh = cafe_shell_likelihood(max);
	int rfsize =  __maxidx(lh,pcafe->rfsize);
	double mlh = lh[rfsize];
   	rfsize += pcafe->rootfamilysizes[0];

	double* pcd = cafe_tree_random_probabilities(cafe_param->pcafe, rfsize, cafe_param->num_random_samples );
	double pv = pvalue( mlh, pcd, cafe_param->num_random_samples);
	pString pstr = cafe_tree_string(pcafe);
	printf("%s\n", pstr->buf );
	printf("Root size: %d with maximum likelihood : %g\n", rfsize, mlh );
	printf("p-value: %g\n", pv);
	memory_free(pcd);
	pcd = NULL;
	string_free(pstr);
}

int cafe_cmd_pvalue(int argc, char* argv[])
{
	int i, j;
	pArrayList pargs = cafe_shell_build_argument(argc, argv);
	pArgument parg;
	FILE* fp;
	if ( ( parg = cafe_shell_get_argument("-o", pargs) ) )
	{
		if ( cafe_param->pcafe == NULL )
		{
			fprintf(stderr, "ERROR(pvalue): You must input newick-formated tree.\n");
			return -1;
		}
		if ( cafe_param->lambda == NULL )
		{
			fprintf(stderr, "ERROR(pvalue): You must specify lambda information.\n");
			return -1;
		}
		if ( cafe_pCD ) arraylist_free( cafe_pCD, free);
		if ( (fp = fopen(parg->argv[0] ,"w") ) == NULL )
		{
			fprintf(stderr, "ERROR(pvalue): Cannot open %s in write mode\n", parg->argv[0] );
			return -1;
		}
		cafe_pCD = cafe_conditional_distribution(cafe_param);
		for ( i = 0 ; i < cafe_pCD->size ; i++ )
		{
			double* data = (double*)cafe_pCD->array[i];
			fprintf(fp, "%10.9lg", data[0] );
			for ( j = 1 ; j < cafe_param->num_random_samples ; j++ )
			{
				fprintf(fp, "\t%10.9lg", data[j] );
			}
			fprintf(fp,"\n");
		}
		fclose(fp);
	}
	else if ( ( parg = cafe_shell_get_argument("-i", pargs) ) )
	{
		if ( cafe_pCD ) arraylist_free( cafe_pCD, free);
		if ( (fp = fopen(parg->argv[0],"r") ) == NULL )
		{
			fprintf(stderr, "ERROR(pvalue): Cannot open %s in read mode\n", parg->argv[0] );
			return -1;
		}
		cafe_log(cafe_param,"Loading p-values ... \n");
		pString pstr = string_new();
		cafe_pCD = arraylist_new(11000);

		while( file_read_line( pstr, fp ) )
		{
			char* buf = pstr->buf;
			double * data = (double*) memory_new(cafe_param->num_random_samples, sizeof(double) );
			int k = 0;
			for ( i = 0, j = 0; buf[i] ; i++ )
			{
				if ( buf[i] == '\t' || buf[i] == '\n' )
				{
					buf[i] = '\0';
					sscanf(&buf[j],"%lg", &data[k++] );
					buf[i] = '\t';
					j = i + 1;
				}
			}
			arraylist_add( cafe_pCD, data );
		}
		string_free(pstr);
		fclose(fp);
		cafe_log(cafe_param,"Done Loading p-values ... \n");
	}
	else if ( ( parg = cafe_shell_get_argument("-idx", pargs) ) )
	{
		if ( cafe_pCD == NULL ) 
		{
			cafe_pCD = cafe_conditional_distribution(cafe_param);
		}
		cafe_shell_set_sizes();

		int idx = 0;

		pCafeTree pcafe = cafe_param->pcafe;
		sscanf( parg->argv[0], "%d", &idx );
		cafe_family_set_size(cafe_param->pfamily,idx,pcafe);
		double* cP = (double*)memory_new( pcafe->rfsize, sizeof(double));
		double* lh = ((pCafeNode)pcafe)->likelihoods;
		cafe_tree_p_values(pcafe, cP, cafe_pCD, cafe_param->num_random_samples);
		for (idx = 0 ; idx < pcafe->rfsize  ; idx++ )
		{
			printf("%d\t%lg\t%lg\n", idx + pcafe->rootfamilysizes[0], lh[idx] , cP[idx] );
		}
		memory_free(cP);
		cP = NULL;
	}
	else if ( pargs->size > 0 )
	{
    STDERR_IF( cafe_param->pcafe == NULL, "You did not specify tree: command 'tree'\n" );
		__cafe_cmd_pvalue_print( cafe_shell_parse_familysize(argc-1,&argv[1] ) );
	}
	else
	{
    STDERR_IF( cafe_param->pcafe == NULL, "You did not specify tree: command 'tree'\n" );
		__cafe_cmd_pvalue_print( cafe_shell_set_familysize() );
	}
	arraylist_free( pargs, free );
	return 0;
}


int cafe_cmd_branchlength(int argc, char* argv[])
{
	STDERR_IF( cafe_param->pcafe == NULL, "ERROR(branchlength): You did not specify tree: command 'tree'\n" );

	pString pstr = cafe_tree_string_with_id( cafe_param->pcafe );
	printf("%s\n", pstr->buf );
	string_free(pstr);
	int err = 0;

	if ( argc == 1 )
	{
		err = cafe_shell_set_branchlength();	
	}
	else if ( argc > 2 )
	{
		err = cafe_shell_parse_branchlength(argc-1,&argv[1]);	
	}
	if ( !err )
	{
		cafe_tree_string_print(cafe_param->pcafe);
	}
	if ( cafe_param->pfamily )
	{
		int i;
		for ( i = 0 ; i < cafe_param->pfamily->flist->size; i++ )
		{
			pCafeFamilyItem pitem = (pCafeFamilyItem)cafe_param->pfamily->flist->array[i];
			pitem->maxlh = -1;
		}
	}
	return err;
}

int cafe_cmd_family(int argc, char* argv[])
{
	STDERR_IF( cafe_param->pfamily == NULL, "ERROR(family): You must load family data first: command 'load'\n");

	int i, idx=0;
	pCafeFamilyItem pitem = NULL;
	pArrayList pargs = cafe_shell_build_argument(argc, argv);
	pArgument parg;

	if ( ( parg = cafe_shell_get_argument("-idx", pargs) ) )
	{
		idx = -1; sscanf( parg->argv[0], "%d", &idx );
	}
	else if ( ( parg = cafe_shell_get_argument("-id", pargs) ) )
	{
		idx = cafe_family_get_index( cafe_param->pfamily, parg->argv[0]);
	}
	else if ( ( parg = cafe_shell_get_argument("-add", pargs) ) )
	{
		pCafeFamily pcf = cafe_param->pfamily;
		if ( pcf == NULL )
		{
			pcf = (pCafeFamily)memory_new( 1, sizeof(CafeFamily));
			pcf->flist = arraylist_new(11000);
			pArrayList nlist = cafe_param->pcafe->super.nlist;
			pcf->num_species = (nlist->size+1)/2;
			pcf->species = (char**) memory_new( pcf->num_species, sizeof(char*));
			pcf->index = (int*) memory_new( pcf->num_species, sizeof(int));
			for ( i = 0 ; i < nlist->size ; i+= 2 )
			{
				pcf->index[i] = i;
				pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
				pcf->species[i] = (char*) memory_new( strlen(pnode->name)+1, sizeof(char) );
				strcpy( pcf->species[i], pnode->name );
			}
			cafe_param->pfamily = pcf;
		}
		pCafeFamilyItem pitem = (pCafeFamilyItem) memory_new(1, sizeof(CafeFamilyItem) );
		pitem->id = (char*) memory_new( strlen( parg->argv[0]) + 1, sizeof(char) );  
		pitem->ref = -1;
		pitem->count = (int*) memory_new( pcf->num_species, sizeof(int) );
		pitem->maxlh = -1;
		pitem->desc = NULL;
		strcpy( pitem->id, parg->argv[0] );
		for ( i = 1 ; i <= pcf->num_species ; i++ )
		{
			sscanf( parg->argv[i], "%d", &pitem->count[pcf->index[i]] );
		}
		arraylist_add( pcf->flist, pitem );
	}
	else if ( ( parg = cafe_shell_get_argument("-del", pargs) ) )
	{
								
	}
	else if ( ( parg = cafe_shell_get_argument("-filter", pargs) ) )
	{
		STDERR_IF( cafe_param->pcafe == NULL, "ERROR(family): You did not specify tree: command 'tree'\n" );
		cafe_family_filter( cafe_param );
		cafe_cmd_print_param(0,NULL);
		return 0;
	}
	if ( idx < 0  )
	{
		fprintf(stderr, "ERROR(family): your request not found\n");
		return -1;
	}
	if ( idx >= cafe_param->pfamily->flist->size )
	{
		fprintf(stderr, "ERROR(family): The index range is from 0 to %d\n", cafe_param->pfamily->flist->size );
		return -1;
	}
	pitem = (pCafeFamilyItem)cafe_param->pfamily->flist->array[idx];
	if ( pitem )
	{
		printf("ID: %s\n", pitem->id );
		printf("Desc: %s\n", pitem->desc );
		for ( i = 0 ; i < cafe_param->pfamily->num_species ; i++ )
		{
			printf("%s: %d\n", cafe_param->pfamily->species[i], pitem->count[i] );
		}
		if ( cafe_param->pcafe && cafe_param->pcafe->pbdc_array ) __cafe_cmd_viterbi_family_print(idx);
	}
	arraylist_free( pargs, free );
	return pitem ? 0 : -1;
}

extern double __cafe_best_lambda_search(double* plambda, void* args);
extern double __cafe_best_lambda_mu_search(double* pparams, void* args);
extern double __cafe_cluster_lambda_search(double* plambda, void* args);
extern double __cafe_cluster_lambda_mu_search(double* pparams, void* args);

double cafe_shell_score()
{
	int i=0;
	double score = 0;
	if (cafe_param->parameterized_k_value > 0) {
		if (cafe_param->num_mus > 0) {
			score = -__cafe_cluster_lambda_mu_search(cafe_param->parameters, (void*)cafe_param);
			// print
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			for( i=0; i<cafe_param->num_lambdas; i++) {
				string_pchar_join_double(buf,",", cafe_param->parameterized_k_value, &cafe_param->parameters[i*cafe_param->parameterized_k_value] );
				cafe_log(cafe_param,"Lambda branch %d: %s\n", i, buf);
				buf[0] = '\0';
			}
			for (i=0; i<cafe_param->num_mus; i++) {
				string_pchar_join_double(buf,",", cafe_param->parameterized_k_value, &cafe_param->parameters[cafe_param->num_lambdas*cafe_param->parameterized_k_value+i*cafe_param->parameterized_k_value]);
				cafe_log(cafe_param,"Mu branch %d: %s \n", i, buf);
				buf[0] = '\0';
			}
			if (cafe_param->parameterized_k_value > 0) {
				string_pchar_join_double(buf,",", cafe_param->parameterized_k_value, cafe_param->k_weights );
				cafe_log(cafe_param, "p : %s\n", buf);
			}
			cafe_log(cafe_param, "Score: %f\n", score);
			
		}
		else {  
			score = -__cafe_cluster_lambda_search(cafe_param->parameters, (void*)cafe_param);
			// print
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			string_pchar_join_double(buf,",", cafe_param->num_lambdas*cafe_param->parameterized_k_value, cafe_param->parameters );
			cafe_log(cafe_param,"Lambda : %s\n", buf);
			buf[0] = '\0';
			if (cafe_param->parameterized_k_value > 0) {
				string_pchar_join_double(buf,",", cafe_param->parameterized_k_value, cafe_param->k_weights );
				cafe_log(cafe_param, "p : %s\n", buf);
			}
			cafe_log(cafe_param, "Score: %f\n", score);		
		}
	}
	else {
		if (cafe_param->num_mus > 0) {
			score = -__cafe_best_lambda_mu_search(cafe_param->parameters, (void*)cafe_param);
			// print
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			string_pchar_join_double(buf,",", cafe_param->num_lambdas, cafe_param->parameters );
			cafe_log(cafe_param,"Lambda : %s ", buf, score);
			buf[0] = '\0';
			string_pchar_join_double(buf,",", cafe_param->num_mus, cafe_param->parameters+cafe_param->num_lambdas );
			cafe_log(cafe_param,"Mu : %s & Score: %f\n", buf, score);		
		}
		else {
			score = -__cafe_best_lambda_search(cafe_param->parameters, (void*)cafe_param);
			// print
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			string_pchar_join_double(buf,",", cafe_param->num_lambdas, cafe_param->parameters );
			cafe_log(cafe_param,"Lambda : %s & Score: %f\n", buf, score);		
		}
	}
	return score;
}

int cafe_cmd_score(int argc, char* argv[])
{
  double score = cafe_shell_score();
  cafe_log(cafe_param,"%lf\n",  score);
  if (cafe_param->parameterized_k_value > 0) {
    cafe_family_print_cluster_membership(cafe_param);
  }
  cafe_shell_set_sizes();
  return 0;
}

int cafe_cmd_retrieve(int argc, char* argv[] )
{
	cafe_shell_clear_param(cafe_param, 1);
	if ( cafe_report_retrieve_data(argv[1], cafe_param)	== -1 )
	{
		return -1;
	}		
	return 0;
}

int __cafe_cmd_log(int argc ,char* argv[] )
{
	if ( cafe_param->str_log )
	{
		string_free( cafe_param->str_log );
		fclose( cafe_param->flog );
		cafe_param->str_log = NULL;
	}
	if ( !strcmp( argv[0], "stdout" ) )
	{
		cafe_param->str_log = NULL;
		cafe_param->flog = stdout;
	}
	else
	{
		cafe_param->str_log = string_join(" ",argc, argv );
		if (  ( cafe_param->flog = fopen( cafe_param->str_log->buf, "a" ) ) == NULL )
		{
			fprintf(stderr, "ERROR(log): Cannot open log file: %s\n", cafe_param->str_log->buf );	
			string_free( cafe_param->str_log );
			cafe_param->flog = stdout;
			return -1;
		}
	}
	return 0;
}

int cafe_cmd_log(int argc, char* argv[] )
{
	if ( argc == 1 ) 
	{
		printf( "Log: %s\n", cafe_param->flog == stdout ? "stdout" : cafe_param->str_log->buf );
	}
	else
	{
		return __cafe_cmd_log(argc-1,&argv[1]);
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
		string_fadd( pstr, "label.rt( btex $\\lambda=%f$ ", pcnode->lambda );
		last -= 0.15;
		string_fadd( pstr, "etex, mid[%d] + (0,%fu));\n",  pnode->id, last );
	}
	return last;
}

extern double cafe_tree_mp_remark(pString pstr, pTree ptree, pMetapostConfig pmc, va_list ap1);

int __cafe_cmd_extinct_count_zero()
{
	int n;
	int cnt_zero = 0;
	tree_clear_reg((pTree)cafe_param->pcafe);
	pArrayList nlist = ((pTree)cafe_param->pcafe)->nlist;
	((pTree)cafe_param->pcafe)->root->reg = 1;
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
/*
	for ( n = 0 ; n < nlist->size ; n++ )
	{
		pTreeNode pnode = (pTreeNode)nlist->array[n];
		pCafeNode pcnode = (pCafeNode)nlist->array[n];
		printf("%d(%s): %d %d\n", n, n%2==0 ? pcnode->super.name : "null", pnode->reg, pcnode->familysize );
	}
	printf(  "--- : %d\n", cnt_zero );
	return 0;
*/
	return cnt_zero;
}

int cafe_cmd_extinct_old(int argc, char* argv[] )
{
	int i, j;
	STDERR_IF( cafe_param->pfamily == NULL, "ERROR(extinct): You did not load family: command 'load'\n" );
	STDERR_IF( cafe_param->pcafe == NULL, "ERROR(extinct): You did not specify tree: command 'tree'\n" );
	STDERR_IF( cafe_param->lambda == NULL, "ERROR(extinct): You did not set the parameters: command 'lambda' or 'lambdamu'\n" );

	pCafeTree pcafe = cafe_param->pcafe;
	int familysize = cafe_param->pfamily->flist->size;
	double* roots_extinct = (double*) memory_new(familysize, sizeof(double));
	int* roots_num = (int*) memory_new(pcafe->rfsize+1,sizeof(int));
	int* roots_size = (int*) memory_new(familysize, sizeof(int));
	cafe_log( cafe_param, "Extinction count from data:\n");

	pCafeNode croot = (pCafeNode)pcafe->super.root;

	for ( i = 0 ; i < familysize ; i++ )
	{
		cafe_family_set_size(cafe_param->pfamily,i, pcafe);
		cafe_tree_viterbi(pcafe);
		roots_extinct[i] = __cafe_cmd_extinct_count_zero();
		roots_size[i] = croot->familysize;
		roots_num[croot->familysize]++;
	}

	int maxsize = roots_num[1];
	for ( i = 2 ; i <= pcafe->rfsize ; i++ )
	{
		if ( maxsize < roots_num[i] ) maxsize = roots_num[i];	
	}

	double* data = (double*) memory_new(maxsize,sizeof(double));

	pHistogram phist = histogram_new(NULL,0,0);
	pHistogram phist_accu = histogram_new( NULL, 0, 0);
	unsigned int accu_sum  = 0;
	for ( i = 1 ; i <= pcafe->rfsize; i++ )
	{
		if ( roots_num[i] )
		{
			int k = 0;
			int sum = 0;
			for ( j = 0 ; j < familysize ; j++ )
			{
				if ( roots_size[j] == i ) 
				{
					data[k++] = roots_extinct[j];
					sum += roots_extinct[j];
				}
			}
			histogram_set_sparse_data(phist,data,k);
			cafe_log( cafe_param, "--------------------------------\n");
			cafe_log( cafe_param, "Root Size: %d\n", i );
			histogram_print(phist, cafe_param->flog);
			if ( cafe_param->flog != stdout ) histogram_print(phist, NULL);
			cafe_log( cafe_param, "Sum : %d\n", sum );
			accu_sum += sum;
		}
	}
	cafe_log( cafe_param, "--------------------------------\n");
	cafe_log( cafe_param, "Total\n");
	histogram_set_sparse_data(phist_accu, roots_extinct, familysize );
	histogram_print(phist_accu, cafe_param->flog);
	if ( cafe_param->flog != stdout ) histogram_print(phist_accu, NULL);
	cafe_log( cafe_param, "Sum : %d\n", accu_sum );

	histogram_free(phist);
	histogram_free(phist_accu);
	
	memory_free(data);
	data = NULL;
	memory_free(roots_extinct);
	roots_extinct = NULL;
	memory_free(roots_size);
	roots_size = NULL;
	memory_free(roots_num);
	roots_num = NULL;
	return 0;
}

int cafe_cmd_sim_extinct(int argc, char* argv[])
{
	STDERR_IF( cafe_param->pcafe == NULL, "ERROR(simextinct): You did not specify tree: command 'tree'\n" );
	STDERR_IF( cafe_param->lambda == NULL, "ERROR(simextinct): You did not set the parameters: command 'lambda' or 'lambdamu'\n" );

	pArrayList pargs = cafe_shell_build_argument(argc, argv);
	pArgument parg;
	int range[2] = { 1, cafe_param->rootfamily_sizes[1] };
	int num_trials = 10000;

/*
	if ( args->argc == 0 )
	{
		fprintf(stderr, "Usage: %s -t #[=10000]\n", argv[0] );
		fprintf(stderr, "-t : number of trials, default 10000"\n);
		arraylist_free(pargs, free);
		return -1;
	}
*/
	cafe_log( cafe_param, "Extinction count from Monte Carole:\n");
	if ( (parg = cafe_shell_get_argument( "-r", pargs) ) )
	{
		if ( index(parg->argv[0], ':' ) )	
		{
			sscanf( parg->argv[0], "%d:%d", &range[0], &range[1] );		
		}
		else
		{
			sscanf( parg->argv[0], "%d", &range[1] );
			range[0] = range[1];
		}
	}
	cafe_log( cafe_param, "root range: %d ~ %d\n", range[0], range[1] );

	if ( (parg = cafe_shell_get_argument( "-t", pargs) ) )
	{
		sscanf( parg->argv[0], "%d", &num_trials);	
	}
	cafe_log( cafe_param, "# trials: %d\n", num_trials );

	if ( range[0] > range[1] || range[1] > cafe_param->rootfamily_sizes[1] )
	{
		fprintf(stderr, "ERROR(simextinct): -r : 1 ~ %d\n", cafe_param->rootfamily_sizes[1] );
		arraylist_free(pargs, free);
		return -1;
	}
	arraylist_free(pargs, free);

	int i, r;
	unsigned int accu_sum = 0;
	pHistogram phist_sim = histogram_new(NULL,0,0);
	pHistogram phist_accu = histogram_new(NULL,0,0);
	double* data = (double*) memory_new( num_trials, sizeof(double) );
	for ( r = range[0] ; r <= range[1] ; r++ )
	{
		int cnt_zero = 0;
		for(  i = 0 ; i < num_trials ; i++ )
		{
			cafe_tree_random_familysize(cafe_param->pcafe, r);
			data[i] = __cafe_cmd_extinct_count_zero();
			cnt_zero += data[i];
		}
		cafe_log( cafe_param, "------------------------------------------\n");
		cafe_log( cafe_param, "Root size: %d\n", r );
		histogram_set_sparse_data(phist_sim,data,num_trials);
		histogram_merge( phist_accu, phist_sim );
		histogram_print(phist_sim,cafe_param->flog );
		if ( cafe_param->flog != stdout ) histogram_print(phist_sim,NULL);
		cafe_log( cafe_param, "Sum : %d\n", cnt_zero );
		accu_sum +=  cnt_zero;
	}

	cafe_log( cafe_param, "------------------------------------------\n");
	cafe_log( cafe_param, "Total\n", r );
	histogram_print(phist_accu,cafe_param->flog );
	if ( cafe_param->flog != stdout ) histogram_print(phist_accu,NULL);
	cafe_log( cafe_param, "Sum : %d\n", accu_sum );


	histogram_free(phist_sim);
	histogram_free(phist_accu);
	tree_clear_reg( (pTree)cafe_param->pcafe );
	memory_free(data);
	data = NULL;
	return 0;
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

int cafe_cmd_root_dist(int argc, char* argv[])
{
	pArrayList pargs = cafe_shell_build_argument(argc, argv);
	pArgument parg;
	STDERR_IF( cafe_param->pcafe == NULL, "ERROR(rootdist): You did not specify tree: command 'tree'\n" );
	
	if ( argc < 2 )
	{
		
		STDERR_IF( cafe_param->pfamily == NULL, "ERROR(rootdist): You did not load family: command 'load'\n" );
		STDERR_IF( cafe_param->lambda == NULL, "ERROR(rootdist): You did not set the parameters: command 'lambda' or 'lambdamu'\n" );
		cafe_log( cafe_param, "-----------------------------------------------------------\n");
		cafe_log( cafe_param, "Family information: %s\n", cafe_param->str_fdata->buf );
		cafe_log( cafe_param, "Log: %s\n", cafe_param->flog == stdout ? "stdout" : cafe_param->str_log->buf );
		if( cafe_param->pcafe ) 
		{
			pString pstr = phylogeny_string( (pTree)cafe_param->pcafe, NULL );
			cafe_log( cafe_param, "Tree: %s\n", pstr->buf ); 
			string_free(pstr);
		}
		if ( cafe_param->lambda )
		{
			pString pstr = cafe_tree_string_with_lambda(cafe_param->pcafe);
			cafe_log( cafe_param, "Lambda: %s\n", pstr->buf );
			string_free(pstr);
		}
		cafe_log( cafe_param, "The number of families is %d\n", cafe_param->pfamily->flist->size );
		int i;
		pCafeTree pcafe = cafe_param->pcafe;
		cafe_set_birthdeath_cache_thread(cafe_param->pcafe, cafe_param->parameterized_k_value, cafe_param->family_sizes, cafe_param->rootfamily_sizes);
		for(i=0; i< cafe_param->pfamily->flist->size; i++) 
		{
			cafe_family_set_size(cafe_param->pfamily,i, pcafe);
			cafe_tree_viterbi(pcafe);
			cafe_log( cafe_param, "%d\n", ((pCafeNode)pcafe->super.root)->familysize);
		}
		cafe_free_birthdeath_cache(pcafe);
		cafe_log( cafe_param, "\n");
	}
	else if ((parg = cafe_shell_get_argument("-i", pargs)))
	{
		pString file = string_join(" ",parg->argc, parg->argv );
		FILE* fp = fopen(file->buf,"r");
		char buf[STRING_BUF_SIZE];
		if ( fp == NULL )
		{
			fprintf( stderr, "Cannot open file: %s\n", file->buf );
			return -1;
		}
		if ( fgets(buf,STRING_BUF_SIZE,fp) == NULL )
		{
			fclose(fp);
			fprintf( stderr, "Empty file: %s\n", file->buf );
			return -1;
		}
		int i=0;
		int max_rootsize = 0;
		string_pchar_chomp(buf);
		pArrayList data = string_pchar_split( buf, ' ');
		pArrayList max = string_pchar_split( data->array[data->size-1], ':');
		max_rootsize = atoi((char*)max->array[1]);
		arraylist_free(data,NULL);
		if (cafe_param->root_dist) { memory_free( cafe_param->root_dist); cafe_param->root_dist = NULL;}
		cafe_param->root_dist = (int*) memory_new(max_rootsize+1,sizeof(int));
		
		cafe_param->pcafe->rootfamilysizes[0] = cafe_param->rootfamily_sizes[0] = 1;
		cafe_param->pcafe->rootfamilysizes[1] = cafe_param->rootfamily_sizes[1] = max_rootsize;
		cafe_param->pcafe->familysizes[0] = cafe_param->family_sizes[0] = 0;
		cafe_param->pcafe->familysizes[1] = cafe_param->family_sizes[1] = max_rootsize*2;
		cafe_param->pcafe->rfsize = cafe_param->rootfamily_sizes[1] - cafe_param->rootfamily_sizes[0] + 1;
		
		for(  i = 0 ; fgets(buf,STRING_BUF_SIZE,fp) ; i++ )	
		{
			string_pchar_chomp(buf);
			data = string_pchar_split( buf, ' ');
			cafe_param->root_dist[atoi(data->array[0])] = atoi(data->array[1]);
		}
		arraylist_free(data,NULL);
	}
	arraylist_free(pargs, free);
	
	return 0;
}

int cafe_cmd_generate_random_family(int argc, char* argv[])
{
	if ( argc == 1 )
	{
		fprintf(stderr,"Usage: %s file\n", argv[0] );
		return -1;
	}
	if ( argv[1] ) {    // check directory exists
		char * dirprefix = strdup(argv[1]);
		char * pch = NULL;
		char * prevpch = NULL; 
		char directory[STRING_BUF_SIZE];
		directory[0] = '\0';
		pch = strtok (dirprefix,"/\0");
		while (pch != NULL)
		{
			if (prevpch) {
				strcat(directory, prevpch);
				strcat(directory, "/");
			}
			prevpch = pch;
			pch = strtok (NULL, "/");
		}
		if (strlen(directory) != 0) {
			struct stat st;
			if(stat(directory,&st) != 0) {
				perror(directory);
				fprintf(stderr, "Please create directory %s before running genfamily.\n", directory);
				return -1;
			}
		}
	}
	
	int num_trials = 1;
	pArrayList pargs = cafe_shell_build_argument(argc, argv);
	pArgument parg;
	if ( (parg = cafe_shell_get_argument( "-t", pargs) ) )
	{
		sscanf( parg->argv[0], "%d", &num_trials);	
	}
	arraylist_free(pargs, free);
	pCafeTree pcafe = cafe_param->pcafe;
	int i, j, n;
	int num_families;
	
	STDERR_IF( cafe_param->pcafe == NULL, "ERROR(genfamily): You did not specify tree: command 'tree'\n" );
	if (cafe_param->root_dist == NULL) {
		STDERR_IF( cafe_param->pfamily == NULL, "ERROR(genfamily): You must either load family data or set root size distribution first: command 'load' or 'rootdist'\n");
		cafe_param->root_dist = (int*)memory_new(pcafe->rfsize+1,sizeof(int));
		num_families = cafe_param->pfamily->flist->size;
		pCafeNode croot = (pCafeNode)pcafe->super.root;
		cafe_param->param_set_func(cafe_param,cafe_param->parameters);
		cafe_set_birthdeath_cache_thread(cafe_param->pcafe, cafe_param->parameterized_k_value, cafe_param->family_sizes, cafe_param->rootfamily_sizes);
		printf("Viterbi\n");
		
		for ( i = 0 ; i < num_families; i++ )
		{
			if ( i % 1000 == 0 )
			{
				printf("%d ...\n", i );
			}
			cafe_family_set_size(cafe_param->pfamily,i, pcafe);
			if (cafe_param->parameterized_k_value > 0) {
				cafe_tree_clustered_viterbi(pcafe);
			}
			else {
				cafe_tree_viterbi(pcafe);
			}
			/*		// use EM instead of viterbi to find maximum likelihood root size
			 double* likelihood = cafe_tree_likelihood(pcafe);		// likelihood of the whole tree = multiplication of likelihood of all nodes
			 cafe_param->ML[i] = __max(likelihood,pcafe->rfsize);			// this part find root size condition with maxlikelihood for each family
			 pCafeFamilyItem pitem = (pCafeFamilyItem)cafe_param->pfamily->flist->array[i];
			 if ( pitem->maxlh < 0 )
			 {
			 pitem->maxlh = __maxidx(likelihood,pcafe->rfsize);	
			 }
			 croot->familysize = pcafe->rootfamilysizes[0]+pitem->maxlh;*/
			cafe_param->root_dist[croot->familysize]++;
		}
	}
	else {
		num_families = 0;
		for ( i = 1 ; i <= pcafe->rfsize ; i++ )
		{
			num_families += cafe_param->root_dist[i];
		}					
		STDERR_IF( cafe_param->lambda == NULL, "ERROR(genfamily): You did not specify lambda: command 'lambda'\n" );
		cafe_log( cafe_param, "Using user defined root size distribution for simulation... \n");
		cafe_param->param_set_func(cafe_param,cafe_param->parameters);
		cafe_set_birthdeath_cache_thread(cafe_param->pcafe, cafe_param->parameterized_k_value, cafe_param->family_sizes, cafe_param->rootfamily_sizes);
	}
	
	int t;
	for ( t = 0 ; t < num_trials ; t++ )
	{
		int id = 1;	
		char buf[STRING_BUF_SIZE];
		sprintf(buf,"%s_%d.tab",argv[1], t+1 );
		FILE* fout = fopen( buf, "w" );
		if ( fout == NULL )
		{
			perror(argv[1]);
			return -1;
		}
		sprintf(buf,"%s_%d.truth",argv[1], t+1 );
		FILE* ftruth = fopen( buf, "w" );
		if ( ftruth == NULL )
		{
			perror(argv[1]);
			return -1;
		}
		fprintf(fout, "DESC\tFID" );
		fprintf(ftruth, "DESC\tFID" );
		pArrayList nlist = pcafe->super.nlist;
		for ( i = 0 ; i < nlist->size ; i+=2 )
		{
			pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
			fprintf(fout,"\t%s", pnode->name );
		}
		fprintf(fout,"\n");
		for ( i = 0 ; i < nlist->size ; i++ )
		{
			pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
			if (pnode->name) {
				fprintf(ftruth,"\t%s", pnode->name );
			}
			else {
				fprintf(ftruth,"\t-%d", i );
			}
		}
		fprintf(ftruth,"\n");
		
		int idx = 0;
		pArrayList k_arr = arraylist_new(num_families);
		if (cafe_param->parameterized_k_value > 0) {
			// assign clusters based on proportion (k_weights)
			int k_i = 0;
			for ( k_i = 0; k_i<cafe_param->parameterized_k_value-1; k_i++) {
				for ( j = 0; j<cafe_param->k_weights[i]*num_families; j++) {
					int* newk = memory_new(1, sizeof(int));
					*newk = k_i;
					arraylist_add(k_arr, (void*)newk);
					idx++;
				}
			}
			for (; idx<num_families; idx++) {
				int* newk = memory_new(1, sizeof(int));
				*newk = k_i;
				arraylist_add(k_arr, (void*)newk);
			}
			// shuffle clusters
			arraylist_shuffle(k_arr);
		}
		idx = 0;
		for ( i = 1 ; i <= pcafe->rfsize ; i++ )
		{
			// iterates along root size distribution, and simulates as many families as the frequency of the root size.
			// need to make it pick random k based on k_weights. 
			if ( cafe_param->root_dist[i] ) 
			{
				for ( j = 0 ; j < cafe_param->root_dist[i] ; j++ )
				{
					// cafe_tree_random_familysize uses birthdeath_matrix to calculate the probabilities.
					// if k > 0 point bd to k_bd[] before running cafe_tree_random_familysize to avoid EXC_BAD_ACCESS		
					if (cafe_param->parameterized_k_value > 0) {
						for ( n = 0 ; n < nlist->size ; n++ )
						{
							pCafeNode pcnode = (pCafeNode)nlist->array[n];
							int* pK = (int*)k_arr->array[idx];
							pcnode->birthdeath_matrix = pcnode->k_bd->array[*pK]; // which k_bd to point to is determined by the randomly shuffled k_arr
						}
					}
					// now randomly sample familysize
					cafe_tree_random_familysize(cafe_param->pcafe, i);
					
					// print test leaves
					if (cafe_param->parameterized_k_value > 0) {
						int* pK = (int*)k_arr->array[idx];
						fprintf(fout, "k%d_root%d\t%d", *pK, i, id );
					}
					else {
						fprintf(fout, "root%d\t%d", i, id );
					}
					for ( n = 0 ; n < nlist->size ; n+=2 )
					{
						pCafeNode pnode = (pCafeNode)nlist->array[n];
						fprintf(fout,"\t%d",pnode->familysize );
					}
					fprintf(fout,"\n");

					// print true leaves and inner nodes
					if (cafe_param->parameterized_k_value > 0) {
						int* pK = (int*)k_arr->array[idx];
						fprintf(ftruth, "k%d_root%d\t%d",*pK, i, id );
					}
					else {
						fprintf(ftruth, "root%d\t%d", i, id );
					}
					for ( n = 0 ; n < nlist->size ; n++ )
					{
						pCafeNode pnode = (pCafeNode)nlist->array[n];
						fprintf(ftruth,"\t%d",pnode->familysize );
					}
					fprintf(ftruth,"\n");
					id++;
					idx++;
				}
			}
		}
		
		arraylist_free(k_arr, free);
		fclose(fout);
		fclose(ftruth);
	}
	cafe_free_birthdeath_cache(pcafe);
	//memory_free( cafe_param->root_dist );
	//cafe_param->root_dist = NULL;
	return 0;
}







int cafe_cmd_extinct(int argc, char* argv[])
{
	STDERR_IF( cafe_param->pfamily == NULL, "ERROR(extinct): You did not load family: command 'load'\n" );
	STDERR_IF( cafe_param->pcafe == NULL, "ERROR(extinct): You did not specify tree: command 'tree'\n" );
	STDERR_IF( cafe_param->lambda == NULL, "ERROR(extinct): You did not set the parameters: command 'lambda' or 'lambdamu'\n" );

	int familysize = cafe_param->pfamily->flist->size;
	pCafeTree pcafe = cafe_param->pcafe;
	double* roots_extinct = (double*) memory_new(familysize,sizeof(double));
	int* roots_size = (int*) memory_new(familysize,sizeof(int));
	int* roots_num = (int*) memory_new(pcafe->rfsize+1,sizeof(int));
	double* avg_extinct = (double*) memory_new(pcafe->rfsize+1,sizeof(double));

	pArrayList pargs = cafe_shell_build_argument(argc, argv);
	pArgument parg;
	int num_trials = 1000;

	cafe_log( cafe_param, ">> Data and Simulation for extinction:\n");

	if ( (parg = cafe_shell_get_argument( "-t", pargs) ) )
	{
		sscanf( parg->argv[0], "%d", &num_trials);	
	}
	arraylist_free(pargs, free);

	cafe_log( cafe_param, "# trials: %d\n", num_trials );

	int i,j;
	int total_extinct = 0;

	fprintf(stderr, "Data ...\n");

	pCafeNode croot = (pCafeNode)pcafe->super.root;
	for ( i = 0 ; i < familysize; i++ )
	{
		cafe_family_set_size(cafe_param->pfamily,i, pcafe);
		cafe_tree_viterbi(pcafe);
		roots_size[i] = croot->familysize;
		roots_extinct[i] = __cafe_cmd_extinct_count_zero();
		total_extinct += roots_extinct[i];
		roots_num[roots_size[i]]++;
	}

	pHistogram phist_tmp = histogram_new(NULL,0,0);
	pHistogram* phist_data = (pHistogram*) memory_new(pcafe->rfsize+1,sizeof(pHistogram));
	pHistogram* phist_sim = (pHistogram*) memory_new(pcafe->rfsize+1,sizeof(pHistogram));
	phist_sim[0] = histogram_new(NULL,0,0);
	phist_data[0] = histogram_new(NULL,0,0);
	int maxsize = roots_num[1];
	for ( i = 1 ; i <= pcafe->rfsize ; i++ ) 
	{
		if ( roots_num[i] > maxsize ) maxsize = roots_num[i];
		if ( roots_num[i] )
		{
			phist_data[i] = histogram_new(NULL,0,0);	
			phist_sim[i] = histogram_new(NULL,0,0);	
		}
	}
	histogram_set_sparse_data(phist_data[0],roots_extinct,familysize);

	double* data = (double*) memory_new( maxsize, sizeof(double));
	double* simdata = (double*) memory_new( familysize, sizeof(double));

	cafe_log(cafe_param,"*******************  DATA  *************************\n");

	for ( i = 1 ; i <= pcafe->rfsize; i++ )
	{
		if ( roots_num[i] )
		{
			int k = 0;
			int sum = 0;
			for ( j = 0 ; j < familysize ; j++ )
			{
				if ( roots_size[j] == i ) 
				{
					data[k++] = roots_extinct[j];
					sum += roots_extinct[j];
				}
			}
			histogram_set_sparse_data(phist_data[i],data,k);
			cafe_log(cafe_param,"--------------------------------\n");
			cafe_log(cafe_param,"Root Size: %d\n", i );
			histogram_print(phist_data[i], cafe_param->flog);
			if ( cafe_param->flog != stdout ) histogram_print(phist_data[i], NULL);
			cafe_log(cafe_param,"Extinct: %d\n", sum );
		}
	}

	fprintf(stderr, "Begin simulation...\n");
	cafe_log(cafe_param,"******************* SIMULATION **********************\n");

	pHistogram** phist_sim_n = (pHistogram**)memory_new_2dim(num_trials, pcafe->rfsize+1, sizeof(pHistogram));

	int t;
	for ( t = 0 ; t < num_trials; t++ )
	{
		if ( t % 100 == 0 && t != 0 )
		{
			fprintf(stderr, "\t%d...\n", t );
		}
		int f = 0;
		for ( i = 1 ; i <= pcafe->rfsize ; i++ )
		{
			if ( roots_num[i] ) 
			{
				int k = 0;
				int sum = 0;
				for ( j = 0 ; j < roots_num[i] ; j++ )
				{
					cafe_tree_random_familysize(cafe_param->pcafe, i);
					simdata[f] = __cafe_cmd_extinct_count_zero();
					data[k++] =simdata[f];
					sum += simdata[f];
					f++;
				}
				avg_extinct[i] += sum;
				histogram_set_sparse_data(phist_tmp,data,k);
				histogram_merge( phist_sim[i], phist_tmp );
				phist_sim_n[t][i] = histogram_new(NULL,0,0);
				histogram_merge( phist_sim_n[t][i], phist_tmp );
			}
		}
		histogram_set_sparse_data(phist_tmp,simdata,familysize);
		histogram_merge( phist_sim[0], phist_tmp );
		phist_sim_n[t][0] = histogram_new(NULL,0,0);
		histogram_merge( phist_sim_n[t][0], phist_tmp);
	}

	double* cnt  = (double*) memory_new(num_trials, sizeof(double));
	double avg_total_extinct = 0;
	for ( i = 1 ; i <= pcafe->rfsize ; i++ )
	{
		if ( phist_sim[i] )
		{
			cafe_log(cafe_param,"--------------------------------\n");
			cafe_log(cafe_param,"Root Size: %d\n", i );
			__hg_print_sim_extinct(phist_sim_n,phist_sim,i,phist_tmp,cnt,num_trials);
			avg_extinct[i] /= num_trials;
			avg_total_extinct += avg_extinct[i];
		}
	}



	cafe_log(cafe_param,"*******************  ALL *************************\n");
	cafe_log(cafe_param,">> DATA\n");
	histogram_print(phist_data[0], cafe_param->flog);
	if ( cafe_param->flog != stdout) histogram_print(phist_data[0], NULL);
	cafe_log(cafe_param, "Total Extinct: %d\n", total_extinct );
	cafe_log(cafe_param,">> SIMULATION\n");
	__hg_print_sim_extinct(phist_sim_n,phist_sim,0,phist_tmp,cnt,num_trials);

	memory_free(cnt);
	cnt = NULL;

	for ( i = 0 ; i < pcafe->rfsize ; i++ )
	{
		if ( phist_data[i] ) 
		{
			histogram_free( phist_data[i] );
			histogram_free( phist_sim[i] );
			for ( t = 0 ; t < num_trials ; t++ )
			{
				histogram_free(phist_sim_n[t][i]);
			}
		}
	}

	histogram_free(phist_tmp);
	memory_free(phist_data);
	phist_data= NULL;
	memory_free(phist_sim);
	phist_sim = NULL;
	memory_free_2dim((void**)phist_sim_n,num_trials, 0, NULL);

	memory_free(roots_size);
	roots_size = NULL;
	memory_free(roots_extinct);
	roots_extinct = NULL;
	memory_free(roots_num);
	roots_num = NULL;
	memory_free(data);
	data  = NULL;
	memory_free(avg_extinct);
	avg_extinct = NULL;
	memory_free(simdata);
	simdata = NULL;
	return 0;
}

int cafe_cmd_lh_test(int argc, char* argv[])
{
	pArrayList pargs = cafe_shell_build_argument(argc, argv);
	pArgument parg;
	double lambda=0;
	char dir[255];
	if ( (parg = cafe_shell_get_argument( "-d", pargs) ) )
	{
		strcpy( dir, parg->argv[0] );
	}
	if ( (parg = cafe_shell_get_argument( "-l", pargs) ) )
	{
		sscanf( parg->argv[0], "%lf", &lambda );	
	}
	char tree[STRING_STEP_SIZE];
	if ( (parg = cafe_shell_get_argument( "-t", pargs) ) )
	{
		strcpy( tree, parg->argv[0] );
	}
	FILE* fout = stdout;
	if ( (parg = cafe_shell_get_argument( "-o", pargs) ) )
	{
		if ( (fout=fopen(parg->argv[0],"w")) == NULL )
		{
			perror( parg->argv[0]);
			return -1;
		}
	}
	arraylist_free(pargs, free);

	pArrayList files = dir_file_list(dir, "tab");
	pString pstr_cafe = phylogeny_string( (pTree)cafe_param->pcafe, NULL );
	int i, j;
	for ( i = 0 ; i < files->size ; i++ )
	{
		char* fname = (char*)files->array[i];
		if ( fname[0] == '.' ) continue;
		cafe_shell_dispatch_commandf("load -i %s/%s -p 0.01 -t 10 -l %s", 
				 dir, fname, cafe_param->str_log ? cafe_param->str_log->buf : "stdout" );
		cafe_shell_dispatch_commandf("tree %s", pstr_cafe->buf );
		cafe_shell_dispatch_commandf("lambda -s -l %lf", lambda );
		fprintf(fout,"\t%lf\t%lf", cafe_get_posterior(cafe_param), cafe_param->lambda[0] );
		cafe_family_reset_maxlh(cafe_param->pfamily);
		cafe_shell_dispatch_commandf("lambda -s -v %lf -t %s", lambda, tree );
		fprintf(fout,"\t%lf", cafe_get_posterior(cafe_param) );
		for ( j = 0 ; j < cafe_param->num_lambdas ; j++ )
		{
			fprintf(fout, "\t%lf", cafe_param->lambda[j] );
		}
		fprintf(fout,"\n");
		fflush(fout);
		cafe_shell_clear_param(cafe_param, 0);
	}
	fclose(fout);
	string_free( pstr_cafe );
	arraylist_free(files, free);	
	return 0;
}




pErrorStruct cafe_shell_create_error_matrix_from_estimate(pErrorMeasure errormeasure);
void cafe_shell_free_errorstruct(pErrorStruct errormodel);


int cafe_shell_read_freq_from_measures(char* file1, char* file2, int* sizeFreq)
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
    int maxFamilySize = cafe_param->family_sizes[1];
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




int cafe_shell_read_error_double_measure(char* error1, char* error2, int** observed_pairs, int maxFamilySize)
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

int cafe_shell_read_error_true_measure(char* errorfile, char* truefile, int** observed_pairs, int maxFamilySize)
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




pErrorMeasure cafe_shell_estimate_error_double_measure(char* error1, char* error2, int b_symmetric, int max_diff, int b_peakzero)
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



pErrorMeasure cafe_shell_estimate_error_true_measure(char* errorfile, char* truefile, int b_symmetric, int max_diff, int b_peakzero)
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




int cafe_shell_write_error_matrix( pErrorStruct errormodel, FILE* fp)
{
    int i, j;
    

    fprintf( fp, "maxcnt:%d\n", errormodel->maxfamilysize );	
    fprintf( fp, "cntdiff");	    
    for (j=errormodel->fromdiff; j<=errormodel->todiff; j++) {
        fprintf( fp," %d", j );	        
    }
    fprintf( fp, "\n");	    
    
    for (j=0; j<=errormodel->maxfamilysize; j++) {
        fprintf( fp, "%d", j);	    
        for (i=errormodel->fromdiff; i<=errormodel->todiff; i++) {
            if (0 <= i+j && i+j <= errormodel->maxfamilysize) { 
                fprintf( fp," %2.2f", errormodel->errormatrix[i+j][j] );	// conditional probability of measuring i+j when true count is j        
            }
            else {
                fprintf( fp," #nan");	        
            }
        }
        fprintf( fp, "\n");
    }
        
    //fclose(fp);

    return 0;
}


// conditional probability of measuring i=familysize when true count is j
int __check_error_model_columnsums(pErrorStruct errormodel) 
{
    int i,j=0;
    int diff = errormodel->todiff;
    
    for (j=0; j<diff; j++) {
        // column j
        double columnsum = 0;
        for (i=0; i<= errormodel->maxfamilysize; i++) { 
            columnsum += errormodel->errormatrix[i][j];
        }
        errormodel->errormatrix[0][j] = errormodel->errormatrix[0][j] + (1-columnsum);
    }
    
    // all other columns
    for (j=diff; j<= errormodel->maxfamilysize-diff; j++) {
        double columnsum = 0;
        for (i=0; i<= errormodel->maxfamilysize; i++) {
            columnsum += errormodel->errormatrix[i][j];
        }
        if (abs(1 - columnsum) > 0.00000000000001) {
            for (i=0; i<= errormodel->maxfamilysize; i++) {
                errormodel->errormatrix[i][j] = errormodel->errormatrix[i][j]/columnsum;
            }
        }
    }
        
    for (j=errormodel->maxfamilysize-diff+1; j<=errormodel->maxfamilysize; j++) {
        // column j
        double columnsum = 0;
        for (i=0; i<= errormodel->maxfamilysize; i++) { 
            columnsum += errormodel->errormatrix[i][j];
        }
        errormodel->errormatrix[errormodel->maxfamilysize][j] = errormodel->errormatrix[errormodel->maxfamilysize][j] + (1-columnsum);

    }
    return 0;
}



pErrorStruct cafe_shell_set_error_matrix(pErrorStruct errormodel, char* speciesname)
{
    int i;
    
    if (cafe_param->pfamily->error_ptr == NULL) {
        cafe_param->pfamily->error_ptr = memory_new(cafe_param->pfamily->num_species, sizeof(pErrorStruct));
    }
    if (speciesname) {
        for (i=0; i<cafe_param->pfamily->num_species; i++) {
            if (string_pchar_cmp_ignore_case(cafe_param->pfamily->species[i], speciesname)) {
                cafe_param->pfamily->error_ptr[i] = errormodel;
                pCafeNode pcnode = (pCafeNode)cafe_param->pcafe->super.nlist->array[cafe_param->pfamily->index[i]];
                pcnode->errormodel = errormodel;
                break;
            }
        }
    }
    else { // '-all' specified instead of speciesname 
        for (i=0; i<cafe_param->pfamily->num_species; i++) {
            cafe_param->pfamily->error_ptr[i] = errormodel;
            pCafeNode pcnode = (pCafeNode)cafe_param->pcafe->super.nlist->array[cafe_param->pfamily->index[i]];
            pcnode->errormodel = errormodel;
        }        
    }
    return errormodel;
    
}

pErrorStruct cafe_shell_create_error_matrix_from_estimate(pErrorMeasure errormeasure)
{
	int i=0;
    int j=0;
    
    // allocate new errormodel
    pErrorStruct errormodel = memory_new(1, sizeof(ErrorStruct));

    errormodel->maxfamilysize = errormeasure->maxFamilySize;
    errormodel->fromdiff = -(errormeasure->model_parameter_diff);
    errormodel->todiff = errormeasure->model_parameter_diff;
    errormodel->errorfilename = NULL;
    errormodel->errormatrix = (double**)memory_new_2dim(errormodel->maxfamilysize+1, errormodel->maxfamilysize+1, sizeof(double));
    
    int total_param_num = 0;
    double* total_params = NULL;
    if (errormeasure->b_symmetric) {
        // symmetric
        double sum = 0;
        total_param_num = errormeasure->model_parameter_number+errormeasure->model_parameter_diff+1;
        total_params = memory_new(total_param_num, sizeof(double));
        total_params[errormeasure->model_parameter_diff] = errormeasure->estimates[0];
        sum = errormeasure->estimates[0];
        for(i=1; i<errormeasure->model_parameter_number; i++) {
            total_params[errormeasure->model_parameter_diff+i] = errormeasure->estimates[i];
            sum += 2*errormeasure->estimates[i];
        } 
        total_params[total_param_num-1] = (1-sum)/(double)((errormeasure->maxFamilySize+1)-(errormeasure->model_parameter_diff*2+1));
        // now fill left side
        for(i=0; i<errormeasure->model_parameter_diff; i++) {
            total_params[i] = total_params[abs(total_param_num-1-1-i)];
        } 
    }   
    else {  
        //asymmetric
        double sum = 0;
        total_param_num = errormeasure->model_parameter_number+1;
        total_params = memory_new(total_param_num, sizeof(double));
        for(i=0; i<errormeasure->model_parameter_number; i++) {
            total_params[i] = errormeasure->estimates[i];
            sum += errormeasure->estimates[i];
        } 
        total_params[total_param_num-1] = (1-sum)/(double)((errormeasure->maxFamilySize+1)-(errormeasure->model_parameter_diff*2+1));
    }

    
    // now fill the error matrix column by column
    for (j=0; j<=errormodel->maxfamilysize; j++) {
        int k = 0;  // k is the index of total_params
        for ( i=0; i<errormodel->fromdiff+j; i++) {
            if (i <= errormodel->maxfamilysize) {
                errormodel->errormatrix[i][j] = total_params[total_param_num-1]; // marginal error probability epsilon
            }
        }      
        for ( i=errormodel->fromdiff+j; i<=errormodel->todiff+j; i++) {
            if (i >= 0 && i <= errormodel->maxfamilysize) {
                errormodel->errormatrix[i][j] = total_params[k]; // conditional probability of measuring i+j when true count is j
            }
            k++;
        }      
        for ( i=errormodel->todiff+j+1; i<=errormodel->maxfamilysize; i++) {
            if (i >= 0) {
                errormodel->errormatrix[i][j] = total_params[total_param_num-1]; // marginal error probability epsilon
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



int cafe_cmd_esterror(int argc, char* argv[]) 
{
	pArrayList pargs = cafe_shell_build_argument(argc, argv);
	pArgument parg;
    //pCafeFamily pcf = cafe_param->pfamily;
    
	if ( argc >= 3 )
	{
        // write errormeasure
        pString outfile = NULL;
        FILE* fp = NULL;

        if ((parg = cafe_shell_get_argument("-o", pargs)))  {
            outfile = string_new_with_string(parg->argv[0]);
            if ( ( fp = fopen( outfile->buf , "w" ) ) == NULL )
            {
                fprintf(stderr,"ERROR(esterror): Cannot open %s in write mode.\n", outfile->buf );  
                return -1;
            }
        }
        else {
            fprintf(stderr,"ERROR(esterror): need to specify the output file \"-o outfile\".\n" );  
            fprintf(stderr,"ERROR(esterror): [-dataerror file1 file2] or [-dataerror file1 -datatrue file2] -o outfile.\n" );  
            return -1;
            
        }
        if ((parg = cafe_shell_get_argument("-dataerror", pargs))) 
        {
            
            pErrorMeasure errormeasure = NULL;
            int b_symmetric = 0;    // default is asymmetric model 
            int max_diff = 2;       // default is a maximum difference of two counts 
            int b_peakzero = 0;     // default is no contraint in shape
            
            if (parg->argc == 2) {
                pString errorfile1 = string_new_with_string(parg->argv[0]);
                pString errorfile2 = string_new_with_string(parg->argv[1]);
                // read model specification 
                if ( (parg = cafe_shell_get_argument( "-symm", pargs) ) )
                {
                    b_symmetric = 1;
                }
                if ( (parg = cafe_shell_get_argument( "-diff", pargs) ) )
                {
                    max_diff = atoi(parg->argv[0]);
                }
                if ( (parg = cafe_shell_get_argument( "-peakzero", pargs) ) )
                {
                    b_peakzero = 1;
                }
                // estimate error matrix                
                errormeasure = cafe_shell_estimate_error_double_measure(errorfile1->buf, errorfile2->buf, b_symmetric, max_diff, b_peakzero);
                string_free(errorfile1);
                string_free(errorfile2);
            }
            else if (parg->argc == 1) {
                pString errorfile1 = string_new_with_string(parg->argv[0]);
                if ((parg = cafe_shell_get_argument("-datatrue", pargs)))
                {
                    pString truefile = string_new_with_string(parg->argv[0]);
                    // read model specification 
                    if ( (parg = cafe_shell_get_argument( "-symm", pargs) ) )
                    {
                        b_symmetric = 1;
                    }
                    if ( (parg = cafe_shell_get_argument( "-diff", pargs) ) )
                    {
                        max_diff = atoi(parg->argv[0]);
                    }
                    if ( (parg = cafe_shell_get_argument( "-peakzero", pargs) ) )
                    {
                        b_peakzero = 1;
                    }
                    // estimate error matrix                                    
                    errormeasure = cafe_shell_estimate_error_true_measure(errorfile1->buf, truefile->buf, b_symmetric, max_diff, b_peakzero);
                    string_free(truefile);
                }
                else {
                        fprintf(stderr,"ERROR(esterror): we need another data file with error or a true data file to compare.\n" );
                        fprintf(stderr,"ERROR(esterror): [-dataerror file1 file2] or [-dataerror file1 -datatrue file2] -o outfile.\n" );
                        fclose(fp);
                        string_free(outfile); 
                        string_free(errorfile1);
                        return -1;
                }
                string_free(errorfile1);
            }
            else
            {
                fprintf(stderr,"ERROR(esterror): [-dataerror file1 file2] or [-dataerror file1 -datatrue file2] -o outfile.\n" );  
                fclose(fp);
                string_free(outfile);                  
                return -1;
            }
            
            // create errormodel based on errormeasure
            pErrorStruct errormodel = cafe_shell_create_error_matrix_from_estimate(errormeasure);
            
            // write errormodel
            cafe_shell_write_error_matrix(errormodel, fp);
            
        }
        fclose(fp);
        string_free(outfile);                  

    }
    else
    {
        fprintf(stderr,"ERROR(esterror): -o outfile [-dataerror file1 file2] or [-dataerror file1 -datatrue file2] -o outfile.\n" );  
        return -1;
    }
    return 0;
    
}


int cafe_shell_set_error_matrix_from_file(char* filename, char* speciesname)
{
	int i=0;
    int j=0;
    char buf[STRING_BUF_SIZE];
    
    // check if error model for filename already exists 
    pErrorStruct errormodel = NULL;
    if (cafe_param->pfamily->errors) {
        for (i=0; i<cafe_param->pfamily->errors->size; i++) {
            pErrorStruct error = (pErrorStruct)cafe_param->pfamily->errors->array[i];
            if ( string_pchar_cmp_ignore_case(error->errorfilename, filename)) {
                errormodel = error;
                break;
            }
        }
    }
    if (errormodel == NULL) 
    {
        // allocate new errormodel
        errormodel = memory_new(1, sizeof(ErrorStruct));
        FILE* fp = fopen(filename,"r");
        if ( fp == NULL )
        {
            fprintf( stderr, "Cannot open file: %s\n", filename );
            return -1;
        }
        if ( fgets(buf,STRING_BUF_SIZE,fp) == NULL )
        {
            fclose(fp);
            fprintf( stderr, "Empty file: %s\n", filename );
            return -1;
        }
		string_pchar_chomp(buf);
		pArrayList data = string_pchar_split( buf, ' ');
		pArrayList max = string_pchar_split( data->array[0], ':');
		errormodel->maxfamilysize = atoi((char*)max->array[1]);
        if (errormodel->maxfamilysize < cafe_param->family_sizes[1]) {
            errormodel->maxfamilysize = cafe_param->family_sizes[1];
        }
		arraylist_free(data,NULL);
        arraylist_free(max, NULL);

        if ( fgets(buf,STRING_BUF_SIZE,fp) != NULL ) {
            string_pchar_chomp(buf);
            pArrayList data = string_pchar_split( buf, ' ');
            errormodel->fromdiff = atoi((char*)data->array[1]);
            errormodel->todiff = atoi((char*)data->array[data->size-1]);
            arraylist_free(data,NULL);        
        }
        errormodel->errorfilename = strdup(filename);       
        errormodel->errormatrix = (double**)memory_new_2dim(errormodel->maxfamilysize+1, errormodel->maxfamilysize+1, sizeof(double));

        i = 0;
        while(fgets(buf,STRING_BUF_SIZE,fp))	
        {
            string_pchar_chomp(buf);
            data = string_pchar_split( buf, ' ');
            if (data->size == (errormodel->todiff-errormodel->fromdiff)+2) {
                while (j && j < atoi(data->array[0])) {
                    // copy previous line's error model for missing lines. 
                    for ( i=errormodel->fromdiff; i<=errormodel->todiff; i++) {
                        
                        if (i+j >= 0 && i+j <= errormodel->maxfamilysize) {
                            errormodel->errormatrix[i+j][j] = errormodel->errormatrix[i+j-1][j-1];
                        }
                    }      
                    i++;
                }
                // read error model and save in matrix row
                int k = 1;  // k is file column index
                for ( i=errormodel->fromdiff; i<=errormodel->todiff; i++) {
                    assert(j == atoi(data->array[0]));
                    if (i+j >= 0 && i+j <= errormodel->maxfamilysize) {
                        errormodel->errormatrix[i+j][j] = atof(data->array[k]);  // conditional probability of measuring i+j when true count is j
                    }
                    k++;
                }      
                j++;
            }
            arraylist_free(data,NULL);
        }
        while (j && j <= errormodel->maxfamilysize) {
            // copy previous line's error model for missing lines till the end of matrix. 
            for ( i=errormodel->fromdiff; i<=errormodel->todiff; i++) {
                if (i+j >= 0 && i+j <= errormodel->maxfamilysize) {
                    errormodel->errormatrix[i+j][j] = errormodel->errormatrix[i+j-1][j-1]; // conditional probability of measuring i+j when true count is j
                }
            }      
            j++;
        }
        
        // now make sure that columns of the error matrix sums to one.
        __check_error_model_columnsums(errormodel);
        
        if (cafe_param->pfamily->errors == NULL) {
            cafe_param->pfamily->errors = arraylist_new(cafe_param->pfamily->num_species);
        }
        arraylist_add(cafe_param->pfamily->errors, errormodel);
        
    }
    if (cafe_param->pfamily->error_ptr == NULL) {
        cafe_param->pfamily->error_ptr = memory_new(cafe_param->pfamily->num_species, sizeof(pErrorStruct));
    }
    if (speciesname) {
        for (i=0; i<cafe_param->pfamily->num_species; i++) {
            if (string_pchar_cmp_ignore_case(cafe_param->pfamily->species[i], speciesname)) {
                cafe_param->pfamily->error_ptr[i] = errormodel;
                pCafeNode pcnode = (pCafeNode)cafe_param->pcafe->super.nlist->array[cafe_param->pfamily->index[i]];
                pcnode->errormodel = errormodel;
                break;
            }
        }
    }
    else { // '-all' specified instead of speciesname 
        for (i=0; i<cafe_param->pfamily->num_species; i++) {
            cafe_param->pfamily->error_ptr[i] = errormodel;
            pCafeNode pcnode = (pCafeNode)cafe_param->pcafe->super.nlist->array[cafe_param->pfamily->index[i]];
            pcnode->errormodel = errormodel;
        }        
    }
    return 0;
}


int cafe_cmd_error_model(int argc, char* argv[])
{
	pArrayList pargs = cafe_shell_build_argument(argc, argv);
	pArgument parg;
	STDERR_IF( cafe_param->pcafe == NULL, "ERROR(errormodel): You did not specify tree: command 'tree'\n" );
	STDERR_IF( cafe_param->pfamily == NULL, "ERROR(errormodel): You did not load family: command 'load'\n" );
    
	
	if ( argc >= 3 )
	{
        if ((parg = cafe_shell_get_argument("-model", pargs)))
        {
            pString file = string_new_with_string(parg->argv[0]);
            if ( (parg = cafe_shell_get_argument( "-sp", pargs) ) )
            {
                int j = 0;
                for ( j = 0 ; j < parg->argc; j++ )
                {
                    pString species = string_new_with_string(parg->argv[j]);
                    cafe_shell_set_error_matrix_from_file(file->buf, species->buf);
                    string_free(species);
                }

            }
            else if ( (parg = cafe_shell_get_argument( "-all", pargs) ) )
            {
                cafe_shell_set_error_matrix_from_file(file->buf, NULL);
            }
            fprintf(stderr,"errormodel: %s set.\n", file->buf );  
            fprintf(stderr,"errormodel: Remember that the rows in the errormodel file have to add up to 1 (rows in the errormodel file correspond to columns in the errormatrix).\n");  
            fprintf(stderr,"errormodel: The program does not check, only renormalizes.\n");  
            string_free(file);
        }

        if (!cafe_param->pfamily->error_ptr || !cafe_param->pfamily->errors) {
            fprintf(stderr, "ERROR(errormodel): we need an error model specified (-model) or two data files.\n");
            return -1;
        }
	}

	arraylist_free(pargs, free);
	
	return 0;
}


int cafe_shell_rm_error_model(char* speciesname)
{
	int i=0;
    if (cafe_param->pfamily->errors) {
        for (i=0; i<cafe_param->pfamily->num_species; i++) {
            if (string_pchar_cmp_ignore_case(cafe_param->pfamily->species[i], speciesname)) {
                cafe_param->pfamily->error_ptr[i] = NULL;
                pCafeNode pcnode = (pCafeNode)cafe_param->pcafe->super.nlist->array[cafe_param->pfamily->index[i]];
                pcnode->errormodel = NULL;
                break;
            }
        }
    }
    return 0;
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

void cafe_shell_free_error_model()
{
    int i;
    for (i=0; i<cafe_param->pfamily->num_species; i++) {
        cafe_shell_rm_error_model(cafe_param->pfamily->species[i]);
    }
    if (cafe_param->pfamily->errors) {
    arraylist_free(cafe_param->pfamily->errors, (freefunc) cafe_shell_free_errorstruct);
    cafe_param->pfamily->errors = NULL;
    }
    if (cafe_param->pfamily->error_ptr) {
    memory_free(cafe_param->pfamily->error_ptr);
    cafe_param->pfamily->error_ptr = NULL;
    }
}

int cafe_cmd_no_error_model(int argc, char* argv[])
{
	pArrayList pargs = cafe_shell_build_argument(argc, argv);
	pArgument parg;
	STDERR_IF( cafe_param->pcafe == NULL, "ERROR(errormodel): You did not specify tree: command 'tree'\n" );
	STDERR_IF( cafe_param->pfamily == NULL, "ERROR(errormodel): You did not load family: command 'load'\n" );
	
	if (( argc >= 3 ) && (parg = cafe_shell_get_argument("-sp", pargs)))
	{
        pString species = string_join("", parg->argc, parg->argv);
        cafe_shell_rm_error_model(species->buf);
	}
    else if (( argc == 2) && (parg = cafe_shell_get_argument("-all", pargs))) 
    {
        cafe_shell_free_error_model();        
    }
	arraylist_free(pargs, free);
	return 0;
}




int __backup_original_count()
{
    int idx = 0;
    pCafeFamily pcf = cafe_param->pfamily;

    if (pcf->countbackup == NULL) {
        pcf->countbackup = memory_new(pcf->flist->size, sizeof(int*));
    }
    for ( idx = 0 ; idx < pcf->flist->size; idx++ )
    {
        pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[idx];
        pcf->countbackup[idx] = pitem->count;
    }
    return 0;
}


int __restore_original_count()
{
    int idx = 0;
    pCafeFamily pcf = cafe_param->pfamily;
    
    if (((pCafeFamilyItem)pcf->flist->array[0])->count != pcf->countbackup[0]) {
        for ( idx = 0 ; idx < pcf->flist->size ; idx++ )
        {
            pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[idx];
            memory_free(pitem->count);
            pitem->count = pcf->countbackup[idx];
        }
        memory_free(pcf->countbackup);
        pcf->countbackup = NULL;
    }
    return 0;
}


int simulate_misclassification(char* filename)
{
    int idx = 0;
    pCafeFamily pcf = cafe_param->pfamily;
    
    // check if error model is set up.
    if (!pcf->error_ptr || !pcf->errors ) {
        fprintf(stderr,"ERROR(simerror): error model is not set up.\nSet error model using command \"errormodel\" prior to running simerror.\n" );  
        return -1;        
    }
    
    // set up file to write simulated data
    if (!filename) {
        filename = strdup("/dev/null");
    }
    FILE* fp = NULL;
    if ( ( fp = fopen( filename , "w" ) ) == NULL )
    {
        fprintf(stderr,"ERROR(simerror): Cannot open %s in write mode.\n", filename );  
        return -1;
    }
    
    char buf[STRING_STEP_SIZE];
    string_pchar_join(buf,"\t", pcf->num_species, pcf->species );
    fprintf(fp,"Desc\tFamily ID\t%s\n", buf );	
    
    
    for ( idx = 0 ; idx < pcf->flist->size; idx++ )
    {
        int n;
        pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[idx];
        int* originalcount = (int*)memory_new_with_init(pcf->num_species, sizeof(int), pitem->count);
        pitem->count = (int*)memory_new(pcf->num_species, sizeof(int));
        for ( n =  0 ; n < pcf->num_species ; n++ )
        {
            int i;
            int j = originalcount[n];     // column idx
            pErrorStruct errormodel = pcf->error_ptr[n];        // right error model
            double* miss_prob = memory_new(errormodel->maxfamilysize+1, sizeof(double));              // misclassification probability
            
            // copy the conditional misclassification probability given true size (column j of errormatrix) 
            for ( i=0; i<=errormodel->maxfamilysize; i++) {
                miss_prob[i] = errormodel->errormatrix[i][j];
            }
            
            // random sampling based on misclassification probability.
            int r;
            double cumul = 0;
            double rnd = unifrnd();					
            for ( r = 0 ; r <= errormodel->maxfamilysize ; r++ )
            {
                cumul += miss_prob[r];
                if ( rnd <= cumul ) break;
            }
            pitem->count[n] = r;
            if (r > errormodel->maxfamilysize) {
                fprintf(stderr,"ERROR(simerror): something gone wrong.\n" );  
                return -1;
            }

        }		
        fprintf(fp,"%s\t%s\t%d", pitem->desc,pitem->id, pitem->count[0]);	
        for ( n =  1 ; n < pcf->num_species ; n++ )
        {
            fprintf(fp,"\t%d", pitem->count[n]);	
        }
        fprintf(fp,"\n");
        free(originalcount);
    }
       
    fclose(fp);    
    return 0;

}

// to check how simex is working 
// generate data with additional error (Ystar) by applying the misclassification matrix to the true data
// ( random sampling of true factors=f, with probabilities following the column=f of the misclassification matrix  )
// then run simex on the variable with error (Ystar), and compare it with the true model (based on Y) and the naive model(based on Ystar). 


// first find the naive estimator by assuming there is no error in the data.
// then we will add even more error and store the estimates for each dataset with increasing error (k = 0.5, 1, 1.5, 2)
// so in the end we have (1+number of lambda) * (number of parameters) estimates. 
// for each i = k we run the estimate j = B number of times. 
// each time add error by applying the misclassification matrix = errormatrix^k
// update the estimates based on error added data, 
// store the mean estimates for all B runs. 
// then predict the estimates at k=-1 
// 
double _cafe_simerror(char* prefix, int repeat) 
{
	int run;
    for (run=0; run<repeat; run++) 
    {
        char trainfile[STRING_BUF_SIZE];
        sprintf(trainfile, "%s_%d.erred", prefix, run);

        __backup_original_count();
        if (simulate_misclassification(trainfile) < 0) {
            fprintf(stderr,"ERROR(simerror): failed simulating misclassification error.\n");  
            return -1;            
        }
        
        // search parameters
/*        if (cafe_param->num_mus > 0) {
            cafe_best_lambda_mu_by_fminsearch(cafe_param, cafe_param->num_lambdas, cafe_param->num_mus, cafe_param->parameterized_k_value);
        }
        else {
            cafe_best_lambda_by_fminsearch(cafe_param, cafe_param->num_lambdas, cafe_param->parameterized_k_value);
        }*/
        __restore_original_count();

    }
    return 0;
}





int cafe_cmd_simerror(int argc, char* argv[])
{
	pArrayList pargs = cafe_shell_build_argument(argc, argv);
    pCafeFamily pcf = cafe_param->pfamily;

	STDERR_IF( cafe_param->pcafe == NULL, "ERROR(simerror): You did not specify tree: command 'tree'\n" );
	STDERR_IF( cafe_param->pfamily == NULL, "ERROR(simerror): You did not load family: command 'load'\n" );
	STDERR_IF( cafe_param->num_params == 0, "ERROR(simerror): You did not specify birth death model: command 'lambda'\n" );
    STDERR_IF( (!pcf->error_ptr || !pcf->errors ), "ERROR(simerror): error model is not set up.\nSet error model using command \"errormodel\" prior to running simerror.\n");

    int i;
    pString prefix = NULL;
    int repeat = 1; 
	for ( i = 0 ; i < pargs->size ; i++ )
	{
		pArgument parg = (pArgument)pargs->array[i];
		
		// Search for whole family 
        if ( !strcmp( parg->opt, "-pre") )
        {
            prefix = string_join("", parg->argc, parg->argv);
        }
		else if ( !strcmp( parg->opt, "-rep" ) )
        {
            repeat = atoi(parg->argv[0]);
        }	
    }
    _cafe_simerror(prefix->buf, repeat);        
	arraylist_free(pargs, free);
	return 0;
}

