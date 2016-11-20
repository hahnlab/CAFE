#ifndef __CAFE_TREE_H__
#define __CAFE_TREE_H__

#include<family.h>
#include<gmatrix.h>
#include <chooseln_cache.h>

typedef struct
{
	char* id;
	int*  count;
	char* desc;
	int   maxlh;
	int   ref;
	double*  lambda;
	double*	 mu;
	double*	 z_membership;
	int   holder;
}CafeFamilyItem;

typedef CafeFamilyItem* pCafeFamilyItem;


typedef enum 
{
	CAFE_REPORT_TEXT = 0,
	CAFE_REPORT_HTML,
	CAFE_REPORT_PDF,
}enumCafeReport;

/****************************************************************************
 * Cafe Tree
****************************************************************************/

extern double cafe_get_posterior(pCafeParam param);
extern double cafe_set_prior_rfsize_empirical(pCafeParam param);
extern pCafeTree cafe_tree_new(char* sztree, int familysizes[], int rootfamilysizes[], double lambda, double mu);
extern pTreeNode cafe_tree_new_empty_node(pTree pcafe);
extern void cafe_tree_set_parameters(pCafeTree pcafe, int familysizes[], int rootfamilysizes[], double lambda);
extern pCafeTree cafe_tree_copy(pCafeTree psrc);
extern pCafeTree cafe_tree_split(pCafeTree pcafe, int idx );
extern void cafe_tree_free(pCafeTree pcafe);
extern void __cafe_tree_free_node(pTree ptree, pTreeNode ptnode, va_list ap);
extern void cafe_tree_string_name(pString pstr, pPhylogenyNode ptnode);
extern pString cafe_tree_string_with_lambda(pCafeTree pcafe);
extern pString cafe_tree_string_with_id(pCafeTree pcafe);
extern pString cafe_tree_string(pCafeTree pcafe);
extern pString cafe_tree_string_with_familysize_lambda(pCafeTree pcafe);
extern pString cafe_tree_string_with_familysize(pCafeTree pcafe);
extern void cafe_tree_string_print(pCafeTree pcafe);
void __cafe_tree_node_compute_likelihood(pTree ptree, pTreeNode ptnode, va_list ap1);
void compute_internal_node_likelihood(pTree ptree, pTreeNode ptnode, struct chooseln_cache *cache);
void compute_leaf_node_likelihood(pTree ptree, pTreeNode ptnode);
void compute_likelihood(pTree ptree, pTreeNode ptnode, struct chooseln_cache* cache);

extern double* cafe_tree_likelihood(pCafeTree pcafe);
extern void cafe_tree_node_free_clustered_likelihoods (pCafeParam param);
extern double** cafe_tree_clustered_likelihood(pCafeTree pcafe); 
extern void cafe_tree_viterbi(pCafeTree pcafe);
extern void cafe_tree_clustered_viterbi(pCafeTree pcafe, int num_likelihoods);
extern void cafe_tree_viterbi_posterior(pCafeTree pcafe, pCafeParam param);
extern double* cafe_tree_random_probabilities(pCafeTree pcafe, int rootFamilysize, int trials );
extern double* cafe_tree_p_values(pCafeTree pcafe, double* p,  pArrayList pconddist, int cdlen);
extern double** cafe_tree_p_values_of_two_trees(pCafeTree pcafe1, pCafeTree pcafe2,    
		 								   double** p,
		                                   pArrayList pconddist1, pArrayList pconddist2,
										   int cdlen );

extern pCafeParam cafe_copy_parameters(pCafeParam psrc);
extern void cafe_free_copy_parameters(pCafeParam param);

pArrayList cafe_tree_conditional_distribution(pCafeTree pcafe, int range_start, int range_end, int num_trials);
extern double cafe_tree_mp_remark(pString str, pTree ptree, pMetapostConfig pmc, va_list ap);
extern double cafe_tree_mp_annotation(pString str, pTreeNode pnode, pMetapostConfig pmc, va_list ap);
extern pMetapostConfig cafe_tree_get_default_mpconfig(int id, double width, double height );
extern pString cafe_tree_metapost(pCafeTree pcafe, int id, char* title, double width, double height);
extern int cafe_tree_random_familysize(pCafeTree pcafe, int rootFamilysize );


/****************************************************************************
 * Cafe Family
****************************************************************************/
pCafeFamily cafe_family_init(pArrayList data);
void cafe_family_add_item(pCafeFamily pcf, pArrayList data);
extern pCafeFamily cafe_family_new(char* file, int bpatcheck);
extern void cafe_family_item_free(pCafeFamilyItem pitem );
extern void cafe_family_free(pCafeFamily pcf);
extern void cafe_family_set_size(pCafeFamily pcf, int idx, pCafeTree pcafe);
extern int cafe_family_set_species_index(pCafeFamily pcf, pCafeTree pcafe );
extern int cafe_family_get_species_index(pCafeFamily pcf, char* speciesname); 
extern void cafe_family_set_size_for_split(pCafeFamily pcf, int idx, pCafeTree pcafe);
extern void cafe_family_set_size_with_family(pCafeFamily pcf, int idx, pCafeTree pcafe );
extern void cafe_family_set_truesize_with_family(pCafeFamily pcf, int idx, pCafeTree pcafe );
extern void cafe_family_set_size_by_species(char* speciesname, int size, pCafeTree pcafe);
extern int cafe_family_get_index(pCafeFamily pcf, char* szid);
extern pCafeFamilyItem cafe_family_get_family_item(pCafeFamily pcf, char* szid );
extern void cafe_family_set_size_with_family_forced(pCafeFamily pcf, int idx, pCafeTree pcafe);
extern void cafe_family_filter( pCafeParam param );
extern int cafe_family_print_cluster_membership(pCafeParam param);
extern void cafe_family_reset_maxlh(pCafeFamily pcf);
extern void cafe_family_read_validate_species(pCafeParam param, char* file);
extern int cafe_family_split_cvfiles_byfamily(pCafeParam param);
extern void cafe_family_clean_cvfiles_byfamily(pCafeParam param);
extern void cafe_family_split_cvfiles_byspecies(pCafeParam param);
extern void cafe_family_clean_cvfiles_byspecies(pCafeParam param); 
extern void cafe_family_read_query_family(pCafeParam param, char* file);

/****************************************************************************
 * Cafe Main
****************************************************************************/

extern void cafe_log(pCafeParam param, const char* msg, ... );
extern void cafe_set_birthdeath_cache(pCafeParam param);
extern void cafe_set_birthdeath_cache_thread(pCafeTree tree, int k_value, int* family_sizes, int* rootfamily_sizes);
extern double* cafe_best_lambda_by_fminsearch(pCafeParam param, int lambda_len, int k);
extern double* cafe_best_lambda_mu_by_fminsearch(pCafeParam param, int lambda_len, int mu_len, int k );
extern double* cafe_best_lambda_mu_eqbg_by_fminsearch(pCafeParam param, int lambda_len, int mu_len );
extern double* cafe_each_best_lambda_by_fminsearch(pCafeParam param, int lambda_len );
extern void cafe_report(pCafeParam param, int method);
extern pArrayList cafe_conditional_distribution(pCafeParam param);
extern void cafe_lambda_set_default(pCafeParam param, double* lambda);

extern void cafe_free_birthdeath_cache(pCafeTree pcafe);
extern pArrayList cafe_viterbi(pCafeParam param, pArrayList pCD);
extern void cafe_branch_cutting(pCafeParam param);
extern void cafe_likelihood_ratio_test(pCafeParam param);
extern pGMatrix cafe_lambda_distribution(pCafeParam param, int numrange, double** range );

extern int cafe_report_retrieve_data(char* file, pCafeParam param);

void initialize_leaf_likelihoods(double **matrix, int num_rows, int range, int familysize, int num_cols, pErrorStruct errormodel);
void reset_k_likelihoods(pCafeNode pcnode, int k, int num_factors);

#define CAFE_VERSION "3.2"

#endif
