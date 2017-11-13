#include "cafe.h"
#include<stdlib.h>
#include<math.h>
#include<mathfunc.h>
#include <chooseln_cache.h>
#include "time.h"

extern void __phylogeny_free_node(pTree ptree, pTreeNode ptnode, va_list ap1);
extern pBirthDeathCacheArray probability_cache;

pTreeNode cafe_tree_new_empty_node(pTree ptree)
{
	pCafeTree pcafe = (pCafeTree)ptree;
	pCafeNode pcnode = 	(pCafeNode)memory_new(1, sizeof(CafeNode) );
	pcnode->likelihoods = (double*)memory_new(pcafe->size_of_factor,sizeof(double));
	pcnode->viterbi = (int*)memory_new(pcafe->size_of_factor,sizeof(int));
	pcnode->birth_death_probabilities.lambda = pcafe->lambda;
	pcnode->birth_death_probabilities.mu = pcafe->mu;
	pcnode->familysize = -1;
	pcnode->errormodel = NULL;
	phylogeny_clear_node((pPhylogenyNode)pcnode);
	return (pTreeNode)pcnode;		
}

void cafe_tree_set_parameters(pCafeTree pcafe, family_size_range* range, double lambda)
{
	int i;
	pArrayList nlist = pcafe->super.nlist;
	copy_range_to_tree(pcafe, range);

	pcafe->lambda = lambda;

	int rsize = pcafe->rfsize;
	int fsize = range->max - range->min  + 1;
	int max_size = rsize > fsize ? rsize : fsize;
	if ( pcafe->size_of_factor < max_size )
	{
		pcafe->size_of_factor = max_size;
		for ( i = 0;  i < nlist->size ; i++ )
		{
			pCafeNode pcnode = (pCafeNode)nlist->array[i];
			pcnode->likelihoods = (double*)memory_realloc(pcnode->likelihoods,pcafe->size_of_factor,sizeof(double));
			pcnode->viterbi = (int*)memory_realloc(pcnode->viterbi,pcafe->size_of_factor,sizeof(int));
		}
	}
}

pTree __cafe_tree_new(tree_func_node_new new_tree_node_func, int size)
{
	pCafeTree pcafe = (pCafeTree)memory_new(1,sizeof(CafeTree));
	pcafe->size_of_factor = size;
	tree_new_fill((pTree)pcafe, cafe_tree_new_empty_node);
	pcafe->super.size = sizeof(CafeTree);
	return (pTree)pcafe;
}


void cafe_tree_parse_node(pTree ptree, pTreeNode ptnode)
{
	pCafeNode pcnode = (pCafeNode)ptnode;
	char* name = pcnode->super.name;
	if ( name == NULL ) return;
	char* familysize = (char*)strchr(name,'_');
	if ( familysize )
	{
		*familysize++ = '\0';
		pcnode->familysize = atoi(familysize);
	}
}

pCafeTree cafe_tree_new(const char* sztree, family_size_range* range, double lambda, double mu)
{
	int rsize = range->root_max - range->root_min  + 1;
	int fsize = range->max - range->min  + 1;

	assert(strlen(sztree) < STRING_BUF_SIZE);
	char buf[STRING_BUF_SIZE];
	strcpy(buf, sztree);	// needs a writable buffer to work with
	pCafeTree pcafe = (pCafeTree)phylogeny_load_from_string(
			             		buf, 
								__cafe_tree_new, 
			             		cafe_tree_new_empty_node, 
								cafe_tree_parse_node, 
								rsize > fsize ? rsize : fsize );
	if (pcafe == NULL) {
		return NULL;
	}

	copy_range_to_tree(pcafe, range);

	pcafe->lambda = lambda;
	pcafe->mu = mu;

	((pCafeNode)pcafe->super.root)->birth_death_probabilities.lambda = lambda;
	((pCafeNode)pcafe->super.root)->birth_death_probabilities.mu = mu;
	tree_build_node_list((pTree)pcafe);
	return pcafe;
}

void __cafe_tree_free_node(pTree ptree, pTreeNode ptnode, va_list ap1)
{
  va_list ap;
  if (ap1) va_copy(ap, ap1);
	pCafeNode pcnode = (pCafeNode)ptnode;
	if ( pcnode->likelihoods ) memory_free(pcnode->likelihoods);	
	pcnode->likelihoods = NULL;
	if ( pcnode->viterbi ) memory_free(pcnode->viterbi);
	pcnode->viterbi = NULL;
	__phylogeny_free_node(ptree,ptnode,ap);		
  if (ap1) va_end(ap1);
}


void cafe_tree_free(pCafeTree pcafe)
{
	if ( pcafe->super.nlist )
	{
		int i;
		for ( i = 0 ; i < pcafe->super.nlist->size ; i++ )
		{
			__cafe_tree_free_node((pTree)pcafe, (pTreeNode)pcafe->super.nlist->array[i], NULL );
		}
	}
	else
	{
		tree_traveral_prefix((pTree)pcafe,__cafe_tree_free_node);
	}
	pTree  ptree = (pTree)pcafe;
	if ( *ptree->count == 0 && ptree->data ) vector_free( ((pVector)ptree->data), free );
	tree_free(ptree);
}

/*******************************************************************************
 *	Tree Output
 *******************************************************************************/

void cafe_tree_string_name(pString pstr, pPhylogenyNode ptnode)
{
	char buf[STRING_BUF_SIZE];
	int familysize =  ((pCafeNode)ptnode)->familysize;
	int idx = 0;
	if ( ptnode->name || familysize >= 0 )
	{
		if ( ptnode->name ) idx = sprintf(buf,"%s", ptnode->name);
		if ( familysize >= 0) sprintf(&buf[idx],"_%d", familysize );
		string_add(pstr,buf);
	}
}

void cafe_tree_string_familysize_lambda(pString pstr, pPhylogenyNode ptnode)
{
	int familysize =  ((pCafeNode)ptnode)->familysize;
	if ( ptnode->branchlength <= 0 ) return;
	if ( ptnode->name ) string_fadd(pstr,"%s", ptnode->name );
	if ( familysize >= 0) string_fadd(pstr,"<%d>", familysize );
	double lambda = ((pCafeNode)ptnode)->birth_death_probabilities.lambda;
	string_fadd(pstr,"_%lf", lambda );
}

void cafe_tree_string_familysize(pString pstr, pPhylogenyNode ptnode)
{
	int familysize =  ((pCafeNode)ptnode)->familysize;
	if ( ptnode->branchlength <= 0 ) return;
	if ( ptnode->name ) string_fadd(pstr,"%s", ptnode->name );
	if ( familysize >= 0) string_fadd(pstr,"< %d >", familysize );
}

void cafe_tree_string_lambda(pString pstr, pPhylogenyNode ptnode)
{
	if ( ptnode->branchlength <= 0 ) return;
	if ( ptnode->name ) string_fadd(pstr,"%s", ptnode->name );
	string_fadd(pstr,"_%lf", ((pCafeNode)ptnode)->birth_death_probabilities.lambda );
}

void cafe_tree_string_id(pString pstr, pPhylogenyNode pnode)
{
	if ( pnode->name ) string_fadd(pstr,"%s", pnode->name );
	string_fadd(pstr,"<%d>", pnode->super.id );
}

pString cafe_tree_string_with_familysize_lambda(pCafeTree pcafe)
{
	return phylogeny_string((pTree)pcafe,cafe_tree_string_familysize_lambda);
}

pString cafe_tree_string_with_lambda(pCafeTree pcafe)
{
	return phylogeny_string((pTree)pcafe,cafe_tree_string_lambda);
}

pString cafe_tree_string_with_familysize(pCafeTree pcafe)
{
	return phylogeny_string((pTree)pcafe,cafe_tree_string_familysize);
}

pString cafe_tree_string_with_id(pCafeTree pcafe)
{
	return phylogeny_string((pTree)pcafe,cafe_tree_string_id);
}


pString cafe_tree_string(pCafeTree pcafe)
{
	return phylogeny_string((pTree)pcafe,cafe_tree_string_name);
}

void cafe_tree_string_print(pCafeTree pcafe)
{
	pString pstr = cafe_tree_string(pcafe);	
	printf("%s\n", pstr->buf );
	string_free(pstr);
}


double cafe_tree_mp_remark(pString pstr, pTree ptree, pMetapostConfig pmc, va_list ap1)
{
  va_list ap;
  va_copy(ap, ap1);
	char* title = va_arg(ap,char*);
	string_add( pstr, "\n% annotation\n");
	string_fadd( pstr, "label( btex  etex, (0.1u, %fu));\n", pmc->height + 0.7 ) ;
	string_fadd( pstr, "label( btex %s etex, (0.1u, %fu));\n", title, pmc->height + 0.5 ) ;
  va_end(ap);
	return pmc->height + 0.5;
}

/**
* \brief Initialize matrix values according to the error model or to defaults
* if familysize < 0, sets the first (range) values of each row to 1, ignoring the others
* otherwise if errormodel is NULL, initializes all values of each row to 0 except for the one indexed by familysize, which is 1
* otherwise, sets each row of matrix to the row of errormatrix indexed by familysize
* I doubt this function is doing what was intended
*/
void initialize_leaf_likelihoods_for_viterbi(double **matrix, int num_rows, int range, int familysize, int num_cols, pErrorStruct errormodel)
{
	for (int i = 0; i < range; i++)
	{
		for (int k = 0; k < num_rows; k++) {
			if (familysize < 0) {
				matrix[k][i] = 1;
			}
			else {
				if (errormodel) {
					memset((void*)matrix[k], 0, num_cols*sizeof(double));
					for (int j = 0; j<num_cols; j++) {
						// conditional probability of measuring i=familysize when true count is j
						matrix[k][j] = errormodel->errormatrix[familysize][j];
					}

				}
				else {
					memset((void*)matrix[k], 0, num_cols*sizeof(double));
					matrix[k][familysize] = 1;
				}
			}
		}
	}

}

/**
* \brief Set likelihood to 1 for actual value, 0 otherwise, or copy values from an existing errormodel
* Copies likelihood values from an errormodel if one exists, 
* otherwise sets all likelihoods to 0 except for familysize, which is set to 1
*/
void initialize_leaf_likelihoods(pTree ptree, pTreeNode ptnode)
{
	pCafeTree pcafe = (pCafeTree)ptree;
	pCafeNode pcnode = (pCafeNode)ptnode;

	if (pcnode->errormodel) {
		memset((void*)pcnode->likelihoods, 0, pcafe->size_of_factor*sizeof(double));
		for (int j = 0; j<pcafe->size_of_factor; j++) {
			// conditional probability of measuring i=familysize when true count is j
			pcnode->likelihoods[j] = pcnode->errormodel->errormatrix[pcnode->familysize][j];
		}

	}
	else {
		// number of likelihoods should be set from the tree's size_of_factor, 
		// therefore the familysize must be less than this
		assert(pcnode->familysize >= 0 && pcnode->familysize < pcafe->size_of_factor);
		memset((void*)pcnode->likelihoods, 0, pcafe->size_of_factor*sizeof(double));
		pcnode->likelihoods[pcnode->familysize] = 1;
	}
}


void compute_viterbis(pCafeNode node, int k, double *factors, int rootfamilysize_start, int rootfamilysize_end, int familysize_start, int familysize_end)
{
	struct square_matrix *bd = node->k_bd->array[k];
	for (int s = rootfamilysize_start, i = 0; s <= rootfamilysize_end; s++, i++)
	{
		double tmp = 0;
		for (int c = familysize_start, j = 0; c <= familysize_end; c++, j++)
		{
			double val = square_matrix_get(bd, s, c);
			tmp = val * node->k_likelihoods[k][j];
			if (tmp > factors[i])
			{
				factors[i] = tmp;
				node->viterbi[i] = j;
			}

		}
	}
}

void compute_child_factor(pCafeTree pcafe, pCafeNode child, family_size_range* range, double *factors)
{
    // p(node=c,child|s) = p(node=c|s)p(child|node=c) integrated over all c
    // remember child likelihood[c]'s never sum up to become 1 because they are likelihoods conditioned on c's.
    // incoming nodes to don't sum to 1. outgoing nodes sum to 1

    if (!child->birthdeath_matrix)
        node_set_birthdeath_matrix(child, probability_cache, pcafe->k);

    assert(child->birthdeath_matrix != NULL);
    square_matrix_multiply(child->birthdeath_matrix, child->likelihoods, range->root_min, range->root_max, range->min, range->max, factors);
}

void compute_internal_node_likelihood(pTree ptree, pTreeNode ptnode)
{
    pCafeTree pcafe = (pCafeTree)ptree;
    pCafeNode pcnode = (pCafeNode)ptnode;

    family_size_range range;
    range.min = pcafe->familysizes[0];
    range.max = pcafe->familysizes[1];
    if (tree_is_root(ptree, ptnode))
    {
        range.root_min = pcafe->rootfamilysizes[0];
        range.root_max = pcafe->rootfamilysizes[1];
    }
    else
    {
        range.root_min = pcafe->familysizes[0];
        range.root_max = pcafe->familysizes[1];
    }

    double *left_factor = memory_new(pcafe->size_of_factor, sizeof(double));
    double *right_factor = memory_new(pcafe->size_of_factor, sizeof(double));

    pCafeNode child1 = (pCafeNode)((pTreeNode)pcnode)->children->head->data;
    pCafeNode child2 = (pCafeNode)((pTreeNode)pcnode)->children->tail->data;

#pragma omp parallel
#pragma omp single nowait
    {
#pragma omp task
        compute_child_factor(pcafe, child1, &range, left_factor);
#pragma omp task
        compute_child_factor(pcafe, child2, &range, right_factor);
#pragma omp taskwait

        int size = range.root_max - range.root_min + 1;
        assert(size <= pcafe->size_of_factor);
        for (int i = 0; i < size; i++)
        {
            pcnode->likelihoods[i] = left_factor[i] * right_factor[i];
        }
    }
    memory_free(left_factor);
    memory_free(right_factor);
}

void free_probabilities(struct probabilities *probs)
{
	if (probs->param_lambdas) 
	{ 
		memory_free(probs->param_lambdas); 
		probs->param_lambdas = NULL; 
	}
	if (probs->param_mus) 
	{ 
		memory_free(probs->param_mus); 
		probs->param_mus = NULL; 
	}

}


void compute_node_likelihoods(pTree ptree, pTreeNode ptnode, va_list ap1)
{
	if (tree_is_leaf(ptnode))
	{
		initialize_leaf_likelihoods(ptree, ptnode);
	}
	else
	{
		compute_internal_node_likelihood(ptree, ptnode);
	}
}

void compute_tree_likelihoods(pCafeTree pcafe)
{
	tree_traveral_postfix((pTree)pcafe, compute_node_likelihoods);
}

double* get_likelihoods(const pCafeTree pcafe)
{
	return ((pCafeNode)pcafe->super.root)->likelihoods;

}

/**
* \brief Initialize node with probability values that it may need.
* If multiple lambdas are set, k_bd is set to an arraylist of matrices with probability values
* In this case, values are set up to the value of num_lambdas
* otherwise the value birthdeath_matrix is used
* probability values are drawn from the cache argument, which should hold a variety
* of possible values
*/
void node_set_birthdeath_matrix(pCafeNode pcnode, pBirthDeathCacheArray cache, int num_lambdas)
{
	if (pcnode->super.branchlength <= 0)
		return;

	if (pcnode->birth_death_probabilities.param_lambdas) {
		if (pcnode->birth_death_probabilities.param_mus) {
			if (num_lambdas > 0) {
				for (int k = 0; k<num_lambdas; k++) {
					struct square_matrix* bd = birthdeath_cache_get_matrix(cache, pcnode->super.branchlength, pcnode->birth_death_probabilities.param_lambdas[k], pcnode->birth_death_probabilities.param_mus[k]);
					arraylist_add(pcnode->k_bd, bd);
				}
			}
			else {
				pcnode->birthdeath_matrix = birthdeath_cache_get_matrix(cache, pcnode->super.branchlength, pcnode->birth_death_probabilities.param_lambdas[0], pcnode->birth_death_probabilities.param_mus[0]);
			}
		}
		else {
			if (num_lambdas > 0) {
				for (int k = 0; k<num_lambdas; k++) {
					struct square_matrix* bd = birthdeath_cache_get_matrix(cache, pcnode->super.branchlength, pcnode->birth_death_probabilities.param_lambdas[k], pcnode->birth_death_probabilities.mu);
					arraylist_add(pcnode->k_bd, bd);
				}
			}
			else {
				pcnode->birthdeath_matrix = birthdeath_cache_get_matrix(cache, pcnode->super.branchlength, pcnode->birth_death_probabilities.param_lambdas[0], pcnode->birth_death_probabilities.mu);
			}
		}
	}
	else {
		pcnode->birthdeath_matrix = birthdeath_cache_get_matrix(cache, pcnode->super.branchlength, pcnode->birth_death_probabilities.lambda, pcnode->birth_death_probabilities.mu);
	}

}

void add_key(pArrayList arr, double branchlength, double lambda, double mu)
{
    for (int i = 0; i < arr->size; ++i)
    {
        struct BirthDeathCacheKey *key = (struct BirthDeathCacheKey *)arraylist_get(arr, i);
        if (key->branchlength == branchlength &&
            key->lambda == lambda &&
            key->mu == mu)
            return;
    }
    struct BirthDeathCacheKey *key = malloc(sizeof(struct BirthDeathCacheKey));
    memset(key, 0, sizeof(struct BirthDeathCacheKey));
    key->branchlength = branchlength;
    key->lambda = lambda;
    key->mu = mu;
    arraylist_add(arr, key);
}

void get_keys_from_node(pCafeNode pcnode, pArrayList arr, int num_lambdas)
{
    if (pcnode->super.branchlength <= 0)
        return;

    if (pcnode->birth_death_probabilities.param_lambdas) {
        if (pcnode->birth_death_probabilities.param_mus) {
            if (num_lambdas > 0) {
                for (int k = 0; k<num_lambdas; k++) {
                    add_key(arr, pcnode->super.branchlength, pcnode->birth_death_probabilities.param_lambdas[k], pcnode->birth_death_probabilities.param_mus[k]);
                }
            }
            else {
                add_key(arr, pcnode->super.branchlength, pcnode->birth_death_probabilities.param_lambdas[0], pcnode->birth_death_probabilities.param_mus[0]);
            }
        }
        else {
            if (num_lambdas > 0) {
                for (int k = 0; k<num_lambdas; k++) {
                    add_key(arr, pcnode->super.branchlength, pcnode->birth_death_probabilities.param_lambdas[k], pcnode->birth_death_probabilities.mu);
                }
            }
            else {
                add_key(arr, pcnode->super.branchlength, pcnode->birth_death_probabilities.param_lambdas[0], pcnode->birth_death_probabilities.mu);
            }
        }
    }
    else {
        add_key(arr, pcnode->super.branchlength, pcnode->birth_death_probabilities.lambda, pcnode->birth_death_probabilities.mu);
    }

}

void do_node_set_birthdeath(pTree ptree, pTreeNode ptnode, va_list ap1)
{
	va_list ap;
	va_copy(ap, ap1);
	pBirthDeathCacheArray cache = va_arg(ap, pBirthDeathCacheArray);
	va_end(ap);

	pCafeTree pcafe = (pCafeTree)ptree;
	node_set_birthdeath_matrix((pCafeNode)ptnode, cache, pcafe->k);
}

void gather_keys(pTree ptree, pTreeNode ptnode, va_list ap1)
{
    va_list ap;
    va_copy(ap, ap1);
    pArrayList arr = va_arg(ap, pArrayList);
    va_end(ap);

    pCafeTree pcafe = (pCafeTree)ptree;
    get_keys_from_node((pCafeNode)ptnode, arr, pcafe->k);
}

/**
*	Set each node's birthdeath matrix based on its values of branchlength, lambdas, and mus
**/
void cafe_tree_set_birthdeath(pCafeTree pcafe, int max_family_size)
{
    pArrayList arr = arraylist_new(40);
    tree_traveral_prefix((pTree)pcafe, gather_keys, arr);

    pBirthDeathCacheArray cache = birthdeath_cache_init(max_family_size);

#pragma omp parallel
#pragma omp for
    for (int i = 0; i < arr->size; ++i)
    {
        struct BirthDeathCacheKey* key = (struct BirthDeathCacheKey*)arraylist_get(arr, i);
        struct square_matrix *matrix = compute_birthdeath_rates(key->branchlength, key->lambda, key->mu, max_family_size);
#pragma omp critical
        hash_table_add(cache->table, key, sizeof(struct BirthDeathCacheKey), matrix, sizeof(struct square_matrix*));
    }

	tree_traveral_prefix((pTree)pcafe, do_node_set_birthdeath, cache);

    // free the cache without deleting the matrices
    void** keys = NULL;
    hash_table_get_keys(cache->table, &keys);
    free(keys);
    hash_table_delete(cache->table);
    memory_free(cache);
}

void cafe_tree_node_copy(pTreeNode psrc, pTreeNode pdest)
{
	pCafeNode pcsrc, pcdest;
	pcsrc = (pCafeNode)psrc;
	pcdest = (pCafeNode)pdest;
	phylogeny_node_copy(psrc,pdest);
	pcdest->birth_death_probabilities.lambda = pcsrc->birth_death_probabilities.lambda;
	pcdest->familysize = pcsrc->familysize;
	pcdest->birthdeath_matrix = pcsrc->birthdeath_matrix;
}

void __cafe_tree_copy_new_fill(pCafeTree psrc, pCafeTree pdest )
{
	pdest->size_of_factor = psrc->size_of_factor;
	pdest->familysizes[0] = psrc->familysizes[0];
	pdest->familysizes[1] = psrc->familysizes[1];
	pdest->rootfamilysizes[0] = psrc->rootfamilysizes[0];
	pdest->rootfamilysizes[1] = psrc->rootfamilysizes[1];
	pdest->lambda = psrc->lambda;
	pdest->rfsize = psrc->rfsize;
}

pCafeTree cafe_tree_copy(pCafeTree psrc)
{
	pCafeTree pcafe = (pCafeTree)tree_copy((pTree)psrc, 
			                    cafe_tree_new_empty_node, 
							    cafe_tree_node_copy );
	__cafe_tree_copy_new_fill(psrc,pcafe);
	tree_build_node_list((pTree)pcafe);
	return pcafe;
}

pCafeTree cafe_tree_split(pCafeTree pcafe, int idx )
{
	pCafeTree psub = (pCafeTree)phylogeny_split_tree((pTree)pcafe,idx, __cafe_tree_free_node );
	if ( psub )
	{
		__cafe_tree_copy_new_fill(pcafe,psub);
	}
	cafe_tree_set_birthdeath(pcafe, probability_cache->maxFamilysize);
	cafe_tree_set_birthdeath(psub, probability_cache->maxFamilysize);
	return psub;
}

/*******************************************************************************
 *	Random family size
 *******************************************************************************/

void __cafe_tree_node_random_familysize(pTree ptree, pTreeNode pnode, va_list ap1)
{
	va_list ap;
	va_copy(ap, ap1);
	int *max = va_arg(ap, int*);
	int max_family_size = va_arg(ap, int);
	va_end(ap);

	if ( tree_is_root(ptree,pnode) ) return;

	double rnd = unifrnd();					
	double cumul = 0;
	pCafeNode pcnode = (pCafeNode)pnode;
	int parent_family_size = ((pCafeNode)pnode->parent)->familysize;
	int c = 0;
	for (; c < max_family_size-1; c++ )
	{
		cumul += square_matrix_get(pcnode->birthdeath_matrix, parent_family_size, c);
		if ( cumul >= rnd ) break;
	}
	pcnode->familysize = c;
	if (*max < pcnode->familysize)
	{
		*max = pcnode->familysize;
	}
}

/**
*	Sets the family size of each node to random value between 0 and the tree's pbdc_array maxFamilySize.
**/
int cafe_tree_random_familysize(pCafeTree pcafe, int rootFamilysize, int maxFamilySize)
{
	int max = 0;
	((pCafeNode)pcafe->super.root)->familysize = rootFamilysize;
	tree_traveral_prefix( (pTree)pcafe, __cafe_tree_node_random_familysize, &max, maxFamilySize);
	return max;
}

void initialize_leaf_likelihood_clustered(pTree ptree, pTreeNode ptnode)
{
    int* rootfamilysizes;
    pCafeTree pcafe = (pCafeTree)ptree;
    pCafeNode pcnode = (pCafeNode)ptnode;

    int s, i, j, k;

    if (tree_is_root(ptree, ptnode->parent))
    {
        rootfamilysizes = pcafe->rootfamilysizes;
    }
    else
    {
        rootfamilysizes = pcafe->familysizes;
    }
    for (s = rootfamilysizes[0], i = 0; s <= rootfamilysizes[1]; s++, i++)
    {
        for (k = 0; k < pcafe->k; k++) {
            if (pcnode->familysize < 0) {
                //fprintf(stderr, "family size not set\n");
                pcnode->k_likelihoods[k][i] = 1;
            }
            else {
                if (pcnode->errormodel) {
                    memset((void*)pcnode->k_likelihoods[k], 0, pcafe->size_of_factor * sizeof(double));
                    for (j = 0; j<pcafe->size_of_factor; j++) {
                        // conditional probability of measuring i=familysize when true count is j
                        pcnode->k_likelihoods[k][j] = pcnode->errormodel->errormatrix[pcnode->familysize][j];
                    }
                }
                else {
                    memset((void*)pcnode->k_likelihoods[k], 0, pcafe->size_of_factor * sizeof(double));
                    pcnode->k_likelihoods[k][pcnode->familysize] = 1;
                }
            }
        }
    }
}

void compute_internal_node_likelihood_clustered(pTree ptree, pTreeNode ptnode)
{
    int* familysizes;
    int* rootfamilysizes;
    pCafeTree pcafe = (pCafeTree)ptree;
    pCafeNode pcnode = (pCafeNode)ptnode;
    double *tree_factors[2];
    tree_factors[0] = memory_new(pcafe->size_of_factor, sizeof(double));
    tree_factors[1] = memory_new(pcafe->size_of_factor, sizeof(double));

    int s, i, k;

    double lambda = -1;
    double mu = -1;
    if (tree_is_root(ptree, ptnode))
    {
        rootfamilysizes = pcafe->rootfamilysizes;
        familysizes = pcafe->familysizes;
    }
    else
    {
        rootfamilysizes = familysizes = pcafe->familysizes;
    }
    int idx;
    double *factors[2] = { NULL, NULL };
    pCafeNode child[2] = { (pCafeNode)((pTreeNode)pcnode)->children->head->data,
        (pCafeNode)((pTreeNode)pcnode)->children->tail->data };
    for (k = 0; k < pcafe->k; k++)
    {
        lambda = pcnode->birth_death_probabilities.param_lambdas[k];
        if (pcnode->birth_death_probabilities.param_mus) {
            mu = pcnode->birth_death_probabilities.param_mus[k];
        }

        // for each child
        for (idx = 0; idx < 2; idx++)
        {
            factors[idx] = tree_factors[idx];
            memset(factors[idx], 0, pcafe->size_of_factor * sizeof(double));
            for (s = rootfamilysizes[0], i = 0; s <= rootfamilysizes[1]; s++, i++)
            {
                for (int c = familysizes[0], j = 0; c <= familysizes[1]; c++, j++)
                {
                    factors[idx][i] += birthdeath_likelihood_with_s_c(s, c, child[idx]->super.branchlength, lambda, mu, NULL) * child[idx]->k_likelihoods[k][j];
                }
            }
        }
        int size = rootfamilysizes[1] - rootfamilysizes[0] + 1;
        for (i = 0; i < size; i++)
        {
            pcnode->k_likelihoods[k][i] = factors[0][i] * factors[1][i];
        }
    }
    memory_free(tree_factors[0]);
    memory_free(tree_factors[1]);
}

void compute_node_clustered_likelihood(pTree ptree, pTreeNode ptnode, va_list unused)
{
    pCafeTree pcafe = (pCafeTree)ptree;
    int maxFamilySize = MAX(pcafe->rootfamilysizes[1], pcafe->familysizes[1]);
    if (!chooseln_is_init())
    {
        chooseln_cache_init(maxFamilySize);
    }
    else if (get_chooseln_cache_size() < maxFamilySize)
    {
        chooseln_cache_resize(maxFamilySize);
    }


    if (tree_is_leaf(ptnode))
    {
        initialize_leaf_likelihood_clustered(ptree, ptnode);
    }
    else
    {
        compute_internal_node_likelihood_clustered(ptree, ptnode);
    }
}


void __cafe_tree_node_compute_clustered_likelihood_using_cache(pTree ptree, pTreeNode ptnode, va_list unused)
{
    pCafeTree pcafe = (pCafeTree)ptree;
    pCafeNode pcnode = (pCafeNode)ptnode;
    double *tree_factors[2];
    tree_factors[0] = memory_new(pcafe->size_of_factor, sizeof(double));
    tree_factors[1] = memory_new(pcafe->size_of_factor, sizeof(double));

    int size;
    int s, c, i, j, k;
    int* rootfamilysizes;
    int* familysizes;
    double** bd = NULL;

    int maxFamilySize = MAX(pcafe->rootfamilysizes[1], pcafe->familysizes[1]);
    if (!chooseln_is_init())
    {
        chooseln_cache_init(maxFamilySize);
    }
    else if (get_chooseln_cache_size() < maxFamilySize)
    {
        chooseln_cache_resize(maxFamilySize);
    }


    if (tree_is_leaf(ptnode))
    {
        if (tree_is_root(ptree, ptnode->parent))
        {
            rootfamilysizes = pcafe->rootfamilysizes;
        }
        else
        {
            rootfamilysizes = pcafe->familysizes;
        }
        for (s = rootfamilysizes[0], i = 0; s <= rootfamilysizes[1]; s++, i++)
        {
            //			pcnode->likelihoods[i] = pcnode->bd[s][pcnode->familysize];
            for (k = 0; k < pcafe->k; k++) {
                if (pcnode->familysize < 0) {
                    //fprintf(stderr, "family size not set\n");
                    pcnode->k_likelihoods[k][i] = 1;
                }
                else {
                    if (pcnode->errormodel) {
                        memset((void*)pcnode->k_likelihoods[k], 0, pcafe->size_of_factor * sizeof(double));
                        for (j = 0; j<pcafe->size_of_factor; j++) {
                            // conditional probability of measuring i=familysize when true count is j
                            pcnode->k_likelihoods[k][j] = pcnode->errormodel->errormatrix[pcnode->familysize][j];
                        }
                    }
                    else {
                        //bd = pcnode->k_bd->array[k];
                        //pcnode->k_likelihoods[k][i] = bd[s][pcnode->familysize];
                        memset((void*)pcnode->k_likelihoods[k], 0, pcafe->size_of_factor * sizeof(double));
                        pcnode->k_likelihoods[k][pcnode->familysize] = 1;
                    }
                }
            }
        }
    }
    else
    {
        if (tree_is_root(ptree, ptnode))
        {
            rootfamilysizes = pcafe->rootfamilysizes;
            familysizes = pcafe->familysizes;
        }
        else
        {
            rootfamilysizes = familysizes = pcafe->familysizes;
        }
        int idx;
        double *factors[2] = { NULL, NULL };
        pCafeNode child[2] = { (pCafeNode)((pTreeNode)pcnode)->children->head->data,
            (pCafeNode)((pTreeNode)pcnode)->children->tail->data };
        for (k = 0; k < pcafe->k; k++)
        {
            // for each child
            for (idx = 0; idx < 2; idx++)
            {
                {
                    factors[idx] = tree_factors[idx];
                    memset(factors[idx], 0, pcafe->size_of_factor * sizeof(double));
                    bd = child[idx]->k_bd->array[k];
                    for (s = rootfamilysizes[0], i = 0; s <= rootfamilysizes[1]; s++, i++)
                    {
                        for (c = familysizes[0], j = 0; c <= familysizes[1]; c++, j++)
                        {
                            factors[idx][i] += bd[s][c] * child[idx]->k_likelihoods[k][j];
                        }
                    }
                }
            }
            size = rootfamilysizes[1] - rootfamilysizes[0] + 1;
            for (i = 0; i < size; i++)
            {
                pcnode->k_likelihoods[k][i] = factors[0][i] * factors[1][i];
            }
        }
    }

    memory_free(tree_factors[0]);
    memory_free(tree_factors[1]);
}

void cafe_tree_node_free_clustered_likelihoods(pCafeParam param)
{
    int i;
    pArrayList nlist = param->pcafe->super.nlist;
    pTree tlambda = param->lambda_tree;
    if (tlambda == NULL)
    {
        for (i = 0; i < nlist->size; i++)
        {
            pCafeNode pcnode = (pCafeNode)nlist->array[i];
            free_probabilities(&pcnode->birth_death_probabilities);
            if (pcnode->k_likelihoods) { memory_free(pcnode->k_likelihoods); pcnode->k_likelihoods = NULL; }
            if (pcnode->k_bd) { arraylist_free(pcnode->k_bd, NULL); pcnode->k_bd = NULL; }
        }
    }
}

double** cafe_tree_clustered_likelihood(pCafeTree pcafe)
{
    if (probability_cache)
    {
        tree_traveral_postfix((pTree)pcafe, __cafe_tree_node_compute_clustered_likelihood_using_cache, NULL);
    }
    else
    {
        tree_traveral_postfix((pTree)pcafe, compute_node_clustered_likelihood);
    }
    return ((pCafeNode)pcafe->super.root)->k_likelihoods;
}
