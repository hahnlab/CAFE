#ifndef __FAMILIY_H__
#define __FAMILIY_H__

#include<tree.h>
#include<birthdeath.h>
#include<pthread.h>
#include<gmatrix.h>

#define FAMILYSIZEMAX	1000
typedef struct 
{
	Tree	super;
	int		familysizes[2], rootfamilysizes[2];	
	double	lambda;
	double mu;
	int branch_params_cnt;
	int k;
	double* factors[2];
	int		size_of_factor;
	int 	rfsize;
	pBirthDeathCacheArray pbdc_array;
}CafeTree;
typedef CafeTree* pCafeTree;


typedef struct
{
    char* errorfilename;
    int fromdiff;
    int todiff;
    int maxfamilysize;
    double** errormatrix;
}ErrorStruct;
typedef ErrorStruct* pErrorStruct;       



typedef struct
{
	PhylogenyNode super;
	double  lambda;
	double	mu;
	double* param_lambdas;
	double* param_mus;
	//double* param_weights;
	double** k_likelihoods;
	double* likelihoods;
	int*    viterbi;
	int	familysize;	
	double** bd;
	pArrayList k_bd;
    pErrorStruct errormodel;
}CafeNode;
typedef CafeNode*	pCafeNode;


typedef struct
{
	char** species;
	int   num_species;
	int*  index;
    pErrorStruct* error_ptr;    // array of ErrorStruct pointers in the same order as species. the pointers point to errors[]. 
	int   max_size;
	pArrayList flist;   // family sizes
    pArrayList errors;  // list of actual ErrorStruct instances
    int** countbackup;  // space to store the real counts while simulating error
}CafeFamily;
typedef CafeFamily* pCafeFamily;

typedef struct tagCafeParam CafeParam;
typedef CafeParam* pCafeParam;
typedef void (*param_func)(pCafeParam param, double* parameters);
//typedef void (*lambda_func)(pCafeParam param, double* lambda);
typedef void (*branchlength_func)(pCafeParam param, int* t);

struct tagCafeParam
{
	FILE *fout, *flog;
	pString str_fdata, str_fout, str_log;

	pCafeTree pcafe;
	pCafeFamily pfamily;
	
	// Deprecated
	//double* lambda_scores;
	//double* lambda_values;

	int eqbg;
	int posterior;
	double* ML;
	double* MAP;
	double* prior_rfsize;
	double** prior_rfsize_by_family;

	double* parameters;
	int num_params;
	param_func param_set_func;
	branchlength_func branchlength_update_func;
	double bl_augment;
	//lambda_func lambda_set_func;

	double* lambda;
	pTree lambda_tree;
	int num_lambdas;

	double* mu;
	pTree mu_tree;
	int num_mus;
	
	int k;
	double* k_weights;
	double** p_z_membership;
	int fixcluster0;
	//double* params;

	double* prior_poisson_lambda;
	
	int checkconv;
	//int* branchlengths_sorted;
    int* old_branchlength;
	double max_branch_length;
    double sum_branch_length;
	int num_branches;
	int family_sizes[2];
	int rootfamily_sizes[2];
	int* root_dist;
	char* cv_species_name;
	pArrayList cv_test_species_list;
	pArrayList cv_test_count_list;
	int cv_fold;

	double pvalue;
	
	int  num_threads;
	int  num_random_samples;

	double** viterbiPvalues;
	int** expandRemainDecrease;
	int** viterbiNodeFamilysizes;
	double* maximumPvalues;
	double* averageExpansion;
	double** cutPvalues;
	double** likelihoodRatios;

	int quiet;
};


/****************************************************************************
 * Cafe Main
****************************************************************************/

extern void thread_run_with_arraylist(int numthreads, void* (*run)(void*), pArrayList pal );
// cafe tree
extern void cafe_tree_set_birthdeath(pCafeTree pcafe);
#endif
