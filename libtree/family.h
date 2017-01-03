#ifndef __FAMILIY_H__
#define __FAMILIY_H__

#include<tree.h>
#include<birthdeath.h>

#define FAMILYSIZEMAX	1000
typedef struct 
{
	Tree	super;
	int		familysizes[2], rootfamilysizes[2];	
	double	lambda;
	double mu;
	int branch_params_cnt;
	int k;
	int		size_of_factor;
	int 	rfsize;
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

struct probabilities
{
	double  lambda;
	double	mu;
	double* param_lambdas;
	double* param_mus;
};

void free_probabilities(struct probabilities *probs);

/** Struct that holds information about a node in a CafeTree. It extends the 
	PhylogenyNode structure which in turn extends the TreeNode structure.
*/
typedef struct
{
	PhylogenyNode super;
	double** k_likelihoods;
	double* likelihoods;
	int*    viterbi;
	
	/**
	* Value representing the size of a single gene family for the species or node
	* Temporary value since CAFE holds in memory counts for many gene families for
	* each node
	*/
	int	familysize;	


	struct probabilities birth_death_probabilities;

	/** Matrix of precalculated values, indexed by the root family size
		and the family size
	*/
	struct square_matrix* birthdeath_matrix;
	pArrayList k_bd;
    pErrorStruct errormodel;
}CafeNode;
typedef CafeNode*	pCafeNode;

typedef struct
{
	/** Number of nodes in the tree */
	int num_nodes;

	/** Number of gene families for which to keep data */
	int num_rows;

	/** Matrix of calculated P values for each node in the tree and each gene family  */
	double** viterbiPvalues;

	/** array of three integers expand, remain, and decrease for each node in the tree relative to its parent */
	int** expandRemainDecrease;

	int** viterbiNodeFamilysizes;
	double* maximumPvalues;
	double* averageExpansion;
	double** cutPvalues;
} viterbi_parameters;

/**
* \brief Structure representing a matrix of values of family sizes
*
* Species array contains a list of names of species. Each item in flist
* holds a pCafeFamilyItem which holds a description and family ID, and
* and an array of integer family sizes in the order of the species given
*/
typedef struct
{
	char** species;				///< Names (ID's) of the species loaded into the family
	int   num_species;			///< Total number of species loaded
	int*  index;				///< indices of the species into the matching \ref CafeTree that was loaded by the user
    pErrorStruct* error_ptr;    ///< array of ErrorStruct pointers in the same order as species. the pointers point to errors[]. 
	int   max_size;
	pArrayList flist;   ///< family sizes
    pArrayList errors;  ///< list of actual ErrorStruct instances
    int** countbackup;  ///< space to store the real counts while simulating error
}CafeFamily;
typedef CafeFamily* pCafeFamily;

typedef struct tagCafeParam CafeParam;
typedef CafeParam* pCafeParam;
typedef void (*param_func)(pCafeParam param, double* parameters);
//typedef void (*lambda_func)(pCafeParam param, double* lambda);

/**
* \brief Singleton structure that holds all of the global data that Cafe acts on.
*
* Initialized at program startup by \ref cafe_shell_init
*/
struct tagCafeParam
{
	FILE *fout, *flog;
	pString str_fdata, str_log;

	/// tree information stored when the user calls the "tree" command
	pCafeTree pcafe;		

	/// family information stored when the user calls the "load" command
	pCafeFamily pfamily;	
	
	int eqbg;
	int posterior;

	/// Max Likelihood - Initialized by the "load" command with the number of families in the table
	double* ML;

	/// root size condition with max likelihood for each family	- Initialized by the "load" command with the number of families in the table
	double* MAP;

	/// prior is a poisson distribution on the root size based on leaves' size
	double* prior_rfsize;

	double* parameters;
	int num_params;
	param_func param_set_func;

	double* lambda;
	pTree lambda_tree;
	int num_lambdas;

	double* mu;
	pTree mu_tree;
	int num_mus;
	
	int parameterized_k_value;
	double* k_weights;
	double** p_z_membership;
	int fixcluster0;

	int checkconv;
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

	double pvalue;
	
	int  num_threads;
	int  num_random_samples;

	viterbi_parameters viterbi;
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
