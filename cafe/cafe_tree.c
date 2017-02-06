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
	char* familysize = (char*)index(name,'_');
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

double cafe_tree_mp_annotation(pString pstr, pTreeNode pnode, pMetapostConfig pmc, va_list ap1)
{
  va_list ap;
  va_copy(ap, ap1);
	pCafeNode pcnode = (pCafeNode)pnode;
    string_add( pstr, ";\n");
	if ( pcnode->familysize >= 0 )
	{
		//char* title = va_arg(ap,char*);
		va_arg(ap,char*); // pass first
		pCafeParam param = va_arg(ap,pCafeParam);
		int fid = va_arg(ap,int);
		string_fadd( pstr, "label.urt( btex %d ", pcnode->familysize, pnode->id );
		if ( param && param->viterbi.viterbiPvalues  && !tree_is_leaf( pnode )  )
		{
			int bidx = 2*(pnode->id/2); 
			if(  param->viterbi.viterbiPvalues[bidx][fid] != -1 )
			{
//				string_fadd( pstr,"\\small{(%4.3f\\%%, %4.3f\\%%)} ", 
				string_fadd( pstr,"\\small{(%4.3f, %4.3f)} ", 
						    param->viterbi.viterbiPvalues[bidx][fid],
							param->viterbi.viterbiPvalues[bidx+1][fid] );
			}
		}
		string_fadd( pstr, "etex, p[%d]);\n", pnode->id );
	}
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
  va_end(ap);
	return last;
}

pMetapostConfig cafe_tree_get_default_mpconfig(int id, double width, double height )
{
	pMetapostConfig pmc = (pMetapostConfig)memory_new(1,sizeof(MetapostConfig));
	pmc->id = id;
	pmc->unit = MP_UNIT_IN;
	pmc->dir = MP_DIR_VERTICAL;
	pmc->shape = MP_SHAPE_RECT | MP_SHAPE_MOST_CENTER;
	pmc->fmod = cafe_tree_mp_annotation;
	pmc->fremark = cafe_tree_mp_remark;
	pmc->width = width;
	pmc->height = height;
	return pmc;
}

pString cafe_tree_metapost(pCafeTree pcafe, int id, char* title, double width, double height )
{
	pMetapostConfig pmc = cafe_tree_get_default_mpconfig(id,width,height);
	pString pstr = phylogeny_to_mp( (pTree)pcafe, pmc, title, NULL, 0);
	memory_free(pmc);
	pmc = NULL;
	return pstr;
}

/*******************************************************************************
 *	Viterbi
 *******************************************************************************/
/* this is for finding the ML path 
instead of summing the likelihood across all possible innernodes
choose the maximum likelihood across all possible innernodes */
void __cafe_tree_node_compute_viterbi(pTree ptree, pTreeNode ptnode, va_list ap1 )
{
	pCafeTree pcafe = (pCafeTree)ptree;
	pCafeNode pcnode = (pCafeNode)ptnode;
	//double lambda = pcafe->lambda;

	double *tree_factors[2];
	tree_factors[0] = memory_new(pcafe->size_of_factor, sizeof(double));
	tree_factors[1] = memory_new(pcafe->size_of_factor, sizeof(double));
	int size;
	int s,c,i,j; 
	int* rootfamilysizes;
	int* familysizes;

	if ( tree_is_leaf(ptnode) )
	{
		if ( tree_is_root(ptree, ptnode->parent) )
		{
			rootfamilysizes = pcafe->rootfamilysizes;
		}
		else
		{
			rootfamilysizes = pcafe->familysizes;
		}
		double* factors = tree_factors[0];
		memset( factors, 0, pcafe->size_of_factor*sizeof(double));
		for ( s = rootfamilysizes[0], i = 0; s <= rootfamilysizes[1] ; s++, i++ )
		{
			if (pcnode->familysize < 0) { 
				//fprintf(stderr, "family size not set\n");
				pcnode->likelihoods[i] = 1;					
				double tmp = 0;
				for( c = pcafe->familysizes[0], j = 0 ; c <= pcafe->familysizes[1] ; c++, j++ )
				{
					tmp = square_matrix_get(pcnode->birthdeath_matrix, s, c);
					if ( tmp > factors[i] )
					{
						factors[i] = tmp;
						pcnode->viterbi[i] = j;
					}
				}
			}
			else {
                if (pcnode->errormodel) {
                    memset((void*)pcnode->likelihoods, 0, pcafe->size_of_factor*sizeof(double));
                    //int start = pcnode->familysize+pcnode->errormodel->fromdiff;
                    //int end = pcnode->familysize+pcnode->errormodel->todiff;
                    //if (start<0) {i=0;j=i-start;} else {i=start;j=0;}
                    for( j=0; j<pcafe->size_of_factor; j++) {
                        // conditional probability of measuring i=familysize when true count is j
                        pcnode->likelihoods[j] = pcnode->errormodel->errormatrix[pcnode->familysize][j];
                    }
                }
                else {
                //pcnode->likelihoods[i] = pcnode->bd[s][pcnode->familysize];
                memset((void*)pcnode->likelihoods, 0, pcafe->size_of_factor*sizeof(double));
                pcnode->likelihoods[pcnode->familysize] = 1;	                    
                }
			}
		}
	}
	else
	{
		if ( tree_is_root(ptree,ptnode) )
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
		for ( idx = 0 ; idx < 2 ; idx++ )
		{
			{
				factors[idx] = tree_factors[idx];
				memset( factors[idx], 0, pcafe->size_of_factor*sizeof(double));
				for( s = rootfamilysizes[0], i = 0 ; s <= rootfamilysizes[1] ; s++, i++ )
				{
					double tmp = 0;
					for( c = familysizes[0], j = 0 ; c <= familysizes[1] ; c++, j++ )
					{
						tmp = square_matrix_get(child[idx]->birthdeath_matrix, s, c) * child[idx]->likelihoods[j];
						if ( tmp > factors[idx][i] )
						{
							factors[idx][i] = tmp;
							child[idx]->viterbi[i] = j;
						}
					}
				}
			}
		}
		size = rootfamilysizes[1] - rootfamilysizes[0] + 1;
		for ( i = 0 ; i < size ; i++ )
		{
			pcnode->likelihoods[i] = factors[0][i] * factors[1][i];
		}
//		printf("%s : %d, %d\n", pcnode->super.name, pcnode->super.branchlength, i );
	}

	memory_free(tree_factors[0]);
	memory_free(tree_factors[1]);
}

void __cafe_tree_node_backtrack_viterbi(pTree ptree, pTreeNode ptnode, va_list ap1 )
{
	pCafeTree pcafe = (pCafeTree)ptree;
	pCafeNode pcnode = (pCafeNode)ptnode;
	if ( tree_is_leaf(ptnode) && (pcnode->familysize >= 0)) return;

	if ( ptree->root == ptnode )
	{
		pcnode->familysize = pcafe->rootfamilysizes[0] + __maxidx( pcnode->likelihoods, pcafe->rfsize );
		/* check family size for all nodes first */
		if(pcnode->familysize>10000) {
			fprintf(stderr,"ERROR: FamilySize larger than bd array size Something wrong\n"); 
			exit(-1);
		}
		/* end check family size for all nodes first */

	}
	else
	{
		pCafeNode pcparent = (pCafeNode)ptnode->parent;
		int base = tree_is_root(ptree,(pTreeNode)pcparent) ? pcafe->rootfamilysizes[0] : pcafe->familysizes[0];
		pcnode->familysize = pcnode->viterbi[pcparent->familysize - base] + pcafe->familysizes[0];
		/* check family size for all nodes first */
		if(pcnode->familysize>10000) {
			fprintf(stderr,"ERROR: FamilySize larger than bd array size Something wrong\n"); 
			exit(-1);
		}
		/* end check family size for all nodes first */
	}
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
		assert(pcnode->familysize < pcafe->size_of_factor);
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

void __cafe_tree_node_compute_clustered_viterbi(pTree ptree, pTreeNode ptnode, va_list ap1 )
{
	va_list ap;
	va_copy(ap, ap1);
	int num_likelihoods = va_arg(ap, int);
	va_end(ap);

	pCafeTree pcafe = (pCafeTree)ptree;
	double *tree_factors[2];
	tree_factors[0] = memory_new(pcafe->size_of_factor, sizeof(double));
	tree_factors[1] = memory_new(pcafe->size_of_factor, sizeof(double));

	pCafeNode pcnode = (pCafeNode)ptnode;
	
	int size;
	int i,k; 
	int* rootfamilysizes;
	int* familysizes;
	
	int maxFamilySize =  MAX( pcafe->rootfamilysizes[1], pcafe->familysizes[1]);
	if ( !chooseln_is_init() ) 
	{
		chooseln_cache_init(maxFamilySize);
	}
	else if ( get_chooseln_cache_size() < maxFamilySize ) 
	{
		chooseln_cache_resize(maxFamilySize);
	}
	
	
	if ( tree_is_leaf(ptnode) )
	{
		if ( tree_is_root(ptree, ptnode->parent) )
		{
			rootfamilysizes = pcafe->rootfamilysizes;
		}
		else
		{
			rootfamilysizes = pcafe->familysizes;
		}

		int range = rootfamilysizes[1] - rootfamilysizes[0] + 1;
		initialize_leaf_likelihoods_for_viterbi(pcnode->k_likelihoods, num_likelihoods, range, pcnode->familysize, pcafe->size_of_factor, pcnode->errormodel);
	}
	else
	{
		if ( tree_is_root(ptree,ptnode) )
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
			for ( idx = 0 ; idx < 2 ; idx++ )
			{
				if (k == 0) {
					factors[idx] = tree_factors[idx];
					memset( factors[idx], 0, pcafe->size_of_factor*sizeof(double));
				}
				compute_viterbis(child[idx], k, factors[idx], rootfamilysizes[0], rootfamilysizes[1], familysizes[0], familysizes[1]);
			}
			size = rootfamilysizes[1] - rootfamilysizes[0] + 1;
			for ( i = 0 ; i < size ; i++ )
			{
				pcnode->k_likelihoods[k][i] = factors[0][i] * factors[1][i];
			}
		}
	}
}




void cafe_tree_viterbi(pCafeTree pcafe)
{
	if ( pcafe->super.postfix )
	{
		pArrayList postfix = pcafe->super.postfix;
		pArrayList prefix = pcafe->super.prefix;
		int i = 0;
		for ( i = 0 ; i < postfix->size; i++ )
		{
			__cafe_tree_node_compute_viterbi((pTree)pcafe, (pTreeNode)postfix->array[i], NULL );
		}
		for ( i = 0 ; i < prefix->size; i++ )
		{
			__cafe_tree_node_backtrack_viterbi((pTree)pcafe, (pTreeNode)prefix->array[i], NULL );
		}
	}
	else
	{
		tree_traveral_postfix((pTree)pcafe, __cafe_tree_node_compute_viterbi);
		tree_traveral_prefix((pTree)pcafe, __cafe_tree_node_backtrack_viterbi);
	}

}


void cafe_tree_clustered_viterbi(pCafeTree pcafe, int num_likelihoods)
{
	tree_traveral_postfix((pTree)pcafe, __cafe_tree_node_compute_clustered_viterbi, num_likelihoods);
	tree_traveral_prefix((pTree)pcafe, __cafe_tree_node_backtrack_viterbi);
}


void cafe_tree_viterbi_posterior(pCafeTree pcafe, pCafeParam param)
{
	if ( pcafe->super.postfix )
	{
		pArrayList postfix = pcafe->super.postfix;
		pArrayList prefix = pcafe->super.prefix;
		int i = 0;
		for ( i = 0 ; i < postfix->size; i++ )
		{
			__cafe_tree_node_compute_viterbi((pTree)pcafe, (pTreeNode)postfix->array[i], NULL );
		}
		if (param->posterior) {
			pCafeNode root = (pCafeNode)postfix->array[postfix->size-1];
			if (tree_is_root((pTree)pcafe, (pTreeNode)root)) {
				for(i = 0; i < param->pcafe->rfsize; i++)	// j: root family size
				{
					// likelihood and posterior both starts from 1 instead of 0 
					root->likelihoods[i] = exp(log(root->likelihoods[i])+log(param->prior_rfsize[i]));	//prior_rfsize also starts from 1
				}				
			}
		}
		for ( i = 0 ; i < prefix->size; i++ )
		{
			__cafe_tree_node_backtrack_viterbi((pTree)pcafe, (pTreeNode)prefix->array[i], NULL );
		}
	}
	else
	{
		tree_traveral_postfix((pTree)pcafe, __cafe_tree_node_compute_viterbi);
		tree_traveral_prefix((pTree)pcafe, __cafe_tree_node_backtrack_viterbi);
	}
	
}

void compute_internal_node_likelihood(pTree ptree, pTreeNode ptnode)
{
	pCafeTree pcafe = (pCafeTree)ptree;
	pCafeNode pcnode = (pCafeNode)ptnode;

	int root_start;
	int root_end;
	int family_start;
	int family_end;
	if (tree_is_root(ptree, ptnode))
	{
		root_start = pcafe->rootfamilysizes[0];
		root_end = pcafe->rootfamilysizes[1];
		family_start = pcafe->familysizes[0];
		family_end = pcafe->familysizes[1];
	}
	else
	{
		root_start = pcafe->familysizes[0];
		root_end = pcafe->familysizes[1];
		family_start = pcafe->familysizes[0];
		family_end = pcafe->familysizes[1];
	}

	double *tree_factors[2];
	tree_factors[0] = memory_new(pcafe->size_of_factor, sizeof(double));
	tree_factors[1] = memory_new(pcafe->size_of_factor, sizeof(double));

	int idx;
	int i = 0;
	double *factors[2] = { NULL, NULL };
	pCafeNode child[2] = { (pCafeNode)((pTreeNode)pcnode)->children->head->data,
		(pCafeNode)((pTreeNode)pcnode)->children->tail->data };
	for (idx = 0; idx < 2; idx++)
	{
		{
			factors[idx] = tree_factors[idx];
			memset(factors[idx], 0, pcafe->size_of_factor*sizeof(double));
			for (int s = root_start, i = 0; s <= root_end; s++, i++)
			{
				for (int c = family_start, j = 0; c <= family_end; c++, j++)
				{
					if (!child[idx]->birthdeath_matrix)
						node_set_birthdeath_matrix(child[idx], probability_cache, pcafe->k);

					factors[idx][i] += square_matrix_get(child[idx]->birthdeath_matrix, s, c) * child[idx]->likelihoods[j];		
					// p(node=c,child|s) = p(node=c|s)p(child|node=c) integrated over all c
					// remember child likelihood[c]'s never sum up to become 1 because they are likelihoods conditioned on c's.
					// incoming nodes to don't sum to 1. outgoing nodes sum to 1
				}
			}
		}
	}
	int size = root_end - root_start + 1;
	for (i = 0; i < size; i++)
	{
		pcnode->likelihoods[i] = factors[0][i] * factors[1][i];
	}
	memory_free(tree_factors[0]);
	memory_free(tree_factors[1]);
}

void __cafe_tree_node_compute_clustered_likelihood(pTree ptree, pTreeNode ptnode, va_list unused)
{
	pCafeTree pcafe = (pCafeTree)ptree;
	pCafeNode pcnode = (pCafeNode)ptnode;
	double *tree_factors[2];
	tree_factors[0] = memory_new(pcafe->size_of_factor, sizeof(double));
	tree_factors[1] = memory_new(pcafe->size_of_factor, sizeof(double));

	int size;
	int s,c,i,j,k; 
	int* rootfamilysizes;
	int* familysizes;
	double lambda = -1;
	double mu = -1;
	//double branchlength = pcnode->super.branchlength;	
	
	int maxFamilySize =  MAX( pcafe->rootfamilysizes[1], pcafe->familysizes[1]);
	if ( !chooseln_is_init() ) 
	{
		chooseln_cache_init(maxFamilySize);
	}
	else if ( get_chooseln_cache_size() < maxFamilySize ) 
	{
		chooseln_cache_resize(maxFamilySize);
	}

	
	if ( tree_is_leaf(ptnode) )
	{
		if ( tree_is_root(ptree, ptnode->parent) )
		{
			rootfamilysizes = pcafe->rootfamilysizes;
		}
		else
		{
			rootfamilysizes = pcafe->familysizes;
		}
		for ( s = rootfamilysizes[0], i = 0; s <= rootfamilysizes[1] ; s++, i++ )
		{
			for (k = 0; k < pcafe->k; k++) { 
				lambda = pcnode->birth_death_probabilities.param_lambdas[k];
				if (pcnode->birth_death_probabilities.param_mus) {
					mu = pcnode->birth_death_probabilities.param_mus[k];
				}
				if (pcnode->familysize < 0) { 
					//fprintf(stderr, "family size not set\n");
					pcnode->k_likelihoods[k][i] = 1;					
				}
				else {
                    if (pcnode->errormodel) {
                        memset((void*)pcnode->k_likelihoods[k], 0, pcafe->size_of_factor*sizeof(double));
                        for( j=0; j<pcafe->size_of_factor; j++) {
                            // conditional probability of measuring i=familysize when true count is j
                            pcnode->k_likelihoods[k][j] = pcnode->errormodel->errormatrix[pcnode->familysize][j];
                        }
                    }
                    else {
                    memset((void*)pcnode->k_likelihoods[k], 0, pcafe->size_of_factor*sizeof(double));
                    pcnode->k_likelihoods[k][pcnode->familysize] = 1;	                    
                    }
				}
			}
		}
	}
	else
	{
		if ( tree_is_root(ptree,ptnode) )
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
			for ( idx = 0 ; idx < 2 ; idx++ )
			{
				{	
					factors[idx] = tree_factors[idx];
					memset( factors[idx], 0, pcafe->size_of_factor*sizeof(double));
					for( s = rootfamilysizes[0], i = 0 ; s <= rootfamilysizes[1] ; s++, i++ )
					{
						for( c = familysizes[0], j = 0 ; c <= familysizes[1] ; c++, j++ )
						{
							factors[idx][i] += birthdeath_likelihood_with_s_c(s, c, child[idx]->super.branchlength, lambda, mu, NULL) * child[idx]->k_likelihoods[k][j];
						}
					}
				}
			}
			size = rootfamilysizes[1] - rootfamilysizes[0] + 1;
			for ( i = 0 ; i < size ; i++ )
			{
				pcnode->k_likelihoods[k][i] = factors[0][i] * factors[1][i];
			}
		}
	}
	memory_free(tree_factors[0]);
	memory_free(tree_factors[1]);
}


void __cafe_tree_node_compute_clustered_likelihood_using_cache(pTree ptree, pTreeNode ptnode, va_list unused)
{
	pCafeTree pcafe = (pCafeTree)ptree;
	pCafeNode pcnode = (pCafeNode)ptnode;
	double *tree_factors[2];
	tree_factors[0] = memory_new(pcafe->size_of_factor, sizeof(double));
	tree_factors[1] = memory_new(pcafe->size_of_factor, sizeof(double));

	int size;
	int s,c,i,j,k; 
	int* rootfamilysizes;
	int* familysizes;
	double** bd = NULL;
	
	int maxFamilySize =  MAX( pcafe->rootfamilysizes[1], pcafe->familysizes[1]);
	if ( !chooseln_is_init() ) 
	{
		chooseln_cache_init(maxFamilySize);
	}
	else if ( get_chooseln_cache_size() < maxFamilySize ) 
	{
		chooseln_cache_resize(maxFamilySize);
	}

	
	if ( tree_is_leaf(ptnode) )
	{
		if ( tree_is_root(ptree, ptnode->parent) )
		{
			rootfamilysizes = pcafe->rootfamilysizes;
		}
		else
		{
			rootfamilysizes = pcafe->familysizes;
		}
		for ( s = rootfamilysizes[0], i = 0; s <= rootfamilysizes[1] ; s++, i++ )
		{
//			pcnode->likelihoods[i] = pcnode->bd[s][pcnode->familysize];
			for (k = 0; k < pcafe->k; k++) { 
				if (pcnode->familysize < 0) { 
					//fprintf(stderr, "family size not set\n");
					pcnode->k_likelihoods[k][i] = 1;					
				}
				else {
                    if (pcnode->errormodel) {
                        memset((void*)pcnode->k_likelihoods[k], 0, pcafe->size_of_factor*sizeof(double));
                        for( j=0; j<pcafe->size_of_factor; j++) {
                            // conditional probability of measuring i=familysize when true count is j
                            pcnode->k_likelihoods[k][j] = pcnode->errormodel->errormatrix[pcnode->familysize][j];
                        }
                    }
                    else {
					//bd = pcnode->k_bd->array[k];
					//pcnode->k_likelihoods[k][i] = bd[s][pcnode->familysize];
                    memset((void*)pcnode->k_likelihoods[k], 0, pcafe->size_of_factor*sizeof(double));
                    pcnode->k_likelihoods[k][pcnode->familysize] = 1;	                    
                    }
				}
			}
		}
	}
	else
	{
		if ( tree_is_root(ptree,ptnode) )
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
			for ( idx = 0 ; idx < 2 ; idx++ )
			{
				{	
					factors[idx] = tree_factors[idx];
					memset( factors[idx], 0, pcafe->size_of_factor*sizeof(double));
					bd = child[idx]->k_bd->array[k];
					for( s = rootfamilysizes[0], i = 0 ; s <= rootfamilysizes[1] ; s++, i++ )
					{
						for( c = familysizes[0], j = 0 ; c <= familysizes[1] ; c++, j++ )
						{
							factors[idx][i] += bd[s][c] * child[idx]->k_likelihoods[k][j];
						}
					}
				}
			}
			size = rootfamilysizes[1] - rootfamilysizes[0] + 1;
			for ( i = 0 ; i < size ; i++ )
			{
				pcnode->k_likelihoods[k][i] = factors[0][i] * factors[1][i];
			}
		}
	}

	memory_free(tree_factors[0]);
	memory_free(tree_factors[1]);
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


void cafe_tree_node_free_clustered_likelihoods (pCafeParam param)
{
	int i;
	pArrayList nlist = param->pcafe->super.nlist;
	pTree tlambda = param->lambda_tree;
	if ( tlambda == NULL )
	{
		for ( i = 0 ; i < nlist->size ; i++ )
		{
			pCafeNode pcnode = (pCafeNode)nlist->array[i];
			free_probabilities(&pcnode->birth_death_probabilities);
			if (pcnode->k_likelihoods) { memory_free(pcnode->k_likelihoods); pcnode->k_likelihoods = NULL;}
			if (pcnode->k_bd) { arraylist_free(pcnode->k_bd, NULL); pcnode->k_bd = NULL; }
		}
	}
}



double** cafe_tree_clustered_likelihood(pCafeTree pcafe) 
{

	if (probability_cache)
	{
		if ( pcafe->super.postfix ) {
			int i = 0;
			pArrayList postfix = ((pTree)pcafe)->postfix;
			for ( i = 0 ; i < postfix->size; i++ )
			{
				__cafe_tree_node_compute_clustered_likelihood_using_cache((pTree)pcafe, (pTreeNode)postfix->array[i], NULL);
			}
		}
		else {
			tree_traveral_postfix((pTree)pcafe, __cafe_tree_node_compute_clustered_likelihood_using_cache, NULL);
		}
	}
	else 	
	{
		if ( pcafe->super.postfix ) {
			int i = 0;
			pArrayList postfix = ((pTree)pcafe)->postfix;
			for ( i = 0 ; i < postfix->size; i++ )
			{
				__cafe_tree_node_compute_clustered_likelihood((pTree)pcafe, (pTreeNode)postfix->array[i], NULL);
			}
		}
		else {
			tree_traveral_postfix((pTree)pcafe, __cafe_tree_node_compute_clustered_likelihood);
		}
	}
	return ((pCafeNode)pcafe->super.root)->k_likelihoods;
		
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

void do_node_set_birthdeath(pTree ptree, pTreeNode ptnode, va_list ap1)
{
	pCafeTree pcafe = (pCafeTree)ptree;
	node_set_birthdeath_matrix((pCafeNode)ptnode, probability_cache, pcafe->k);
}


/**
*	Set each node's birthdeath matrix based on its values of branchlength, lambdas, and mus
**/
void cafe_tree_set_birthdeath(pCafeTree pcafe)
{
	tree_traveral_prefix((pTree)pcafe, do_node_set_birthdeath);
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
	cafe_tree_set_birthdeath(pcafe);
	cafe_tree_set_birthdeath(psub);
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
	va_end(ap);

	if ( tree_is_root(ptree,pnode) ) return;

	double rnd = unifrnd();					
	double cumul = 0;
	pCafeNode pcnode = (pCafeNode)pnode;
	int maxFamilysize = probability_cache->maxFamilysize;
	int s = ((pCafeNode)pnode->parent)->familysize;
	int c = 0;
	for (; c < maxFamilysize ; c++ )
	{
		cumul += square_matrix_get(pcnode->birthdeath_matrix, s, c);
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
int cafe_tree_random_familysize(pCafeTree pcafe, int rootFamilysize )
{
	int max = 0;
	((pCafeNode)pcafe->super.root)->familysize = rootFamilysize;
	tree_traveral_prefix( (pTree)pcafe, __cafe_tree_node_random_familysize, &max);
	return max;
}

void cafe_tree_p_values(pCafeTree pcafe,double* pvalues, pArrayList pconddist, int cdlen)
{
	compute_tree_likelihoods(pcafe);
	double* lh = get_likelihoods(pcafe);
	int s;
	for( s = 0 ; s < pcafe->rfsize; s++ )
	{
		pvalues[s] = pvalue( lh[s] , (double*)pconddist->array[s], cdlen);
	}
}

