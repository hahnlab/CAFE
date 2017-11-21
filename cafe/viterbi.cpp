#include <iostream>

#include "viterbi.h"
#include "Globals.h"
#include "pvalue.h"

extern "C" {
#include <family.h>
#include "cafe.h"
#include <pthread.h>
    extern struct chooseln_cache cache;


}


/**************************************************************************
* Viterbi
**************************************************************************/
void viterbi_parameters_init(viterbi_parameters *viterbi, int nnodes, int nrows)
{
	viterbi->num_nodes = nnodes;
	viterbi->num_rows = nrows;
	viterbi->maximumPvalues = (double*)memory_new(nrows, sizeof(double));
	viterbi->averageExpansion.clear();
	viterbi->expandRemainDecrease.clear();
	viterbi->averageExpansion.resize(nnodes);
	viterbi->expandRemainDecrease.resize(nnodes);
}

void viterbi_set_max_pvalue(viterbi_parameters* viterbi, int index, double val)
{
	assert(index < viterbi->num_rows);
	viterbi->maximumPvalues[index] = val;
}


pthread_mutex_t mutex_cafe_viterbi = PTHREAD_MUTEX_INITIALIZER;

void viterbi_sum_probabilities(viterbi_parameters *viterbi, pCafeTree pcafe, pCafeFamilyItem pitem)
{
    pTree ptree = (pTree)pcafe;
    int nnodes = (ptree->nlist->size - 1) / 2;
    for (int j = 0; j < nnodes; j++)
    {
        pCafeNode pcnode = (pCafeNode)ptree->nlist->array[2 * j + 1];
        pCafeNode child[2] = { (pCafeNode)((pTreeNode)pcnode)->children->head->data,
            (pCafeNode)((pTreeNode)pcnode)->children->tail->data };
        for (int k = 0; k < 2; k++)
        {
            double p = square_matrix_get(child[k]->birthdeath_matrix, pcnode->familysize, child[k]->familysize);
            int node_id = 2 * j + k;
            for (int m = 0; m <= pcafe->familysizes[1]; m++)
            {
                auto key = viterbi_parameters::NodeFamilyKey(node_id, pitem);
                if (square_matrix_get(child[k]->birthdeath_matrix, pcnode->familysize, m) == p)
                {
                    viterbi->viterbiPvalues[key] += square_matrix_get(child[k]->birthdeath_matrix, pcnode->familysize, m) / 2.0;
                }
                else if (square_matrix_get(child[k]->birthdeath_matrix, pcnode->familysize, m) < p)
                {
                    viterbi->viterbiPvalues[key] += square_matrix_get(child[k]->birthdeath_matrix, pcnode->familysize, m);
                }
            }
        }
    }
}

void familysize_sanity_check(pTree ptree)
{
  int nnodes = (ptree->nlist->size - 1) / 2;
  /* check family size for all nodes first */
  for (int j = 0; j < nnodes; j++)
  {
    pCafeNode pcnode = (pCafeNode)ptree->nlist->array[j];
    if (pcnode->familysize>10000) {
      fprintf(stderr, "ERROR: Unreasonably large family size (%d). Cannot continue\n", pcnode->familysize);
      exit(-1);
    }
  }
  /* end check family size for all nodes first */
}

void viterbi_section(pCafeFamily pcf, double pvalue, int num_random_samples, viterbi_parameters *viterbi, int i, pCafeTree pcafe, double *cP, std::vector<std::vector<double> >* pCD)
{
    pTree ptree = (pTree)pcafe;

    cafe_family_set_size_with_family_forced(pcf, i, pcafe);

    std::vector<double> p1(pcafe->rfsize);
    cafe_tree_p_values(pcafe, p1, *pCD, num_random_samples);
    viterbi_set_max_pvalue(viterbi, i, __max(&p1[0], pcafe->rfsize));
    cafe_tree_viterbi(pcafe);
    familysize_sanity_check(ptree);

    pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[i];

    viterbi->compute_size_deltas(ptree, pitem);
    if (viterbi->maximumPvalues[i] > pvalue)
    {
        for (int j = 0; j < ptree->nlist->size - 1; j++)
        {
            pCafeNode pnode = (pCafeNode)ptree->nlist->array[j];
            assert(pnode);
            viterbi->viterbiPvalues[viterbi_parameters::NodeFamilyKey(pnode->super.super.id, pitem)] = -1;
        }
        return;
    }

    viterbi_sum_probabilities(viterbi, pcafe, pitem);
}


void* __cafe_viterbi_thread_func(void* ptr)
{
	ViterbiParam* pv = (ViterbiParam *)ptr;

    std::vector<std::vector<double> >* pCD = pv->pCD;
	pCafeTree pcafe = cafe_tree_copy(pv->pcafe);
	int fsize = pv->pfamily->flist->size;
	double* cP = (double*)memory_new(pcafe->rfsize, sizeof(double));
#ifdef VERBOSE
	printf("VITERBI: from %d\n", pv->from);
#endif
	for (int i = pv->from; i < fsize; i += pv->num_threads)
	{
		viterbi_section(pv->pfamily, pv->pvalue, pv->num_random_samples, pv->viterbi, i, pcafe, cP, pCD);
	}
	memory_free(cP);
	cP = NULL;
	cafe_tree_free(pcafe);
	return (NULL);
}

void cafe_viterbi(Globals& globals, viterbi_parameters& viterbi, std::vector<std::vector<double> >* pCD)
{
	pCafeParam param = &globals.param;
	cafe_log(param, "Running Viterbi algorithm....\n");

	ViterbiParam* ptparam = (ViterbiParam*)memory_new(param->num_threads, sizeof(ViterbiParam));
	pTree ptree = (pTree)param->pcafe;

	int nrows = param->pfamily->flist->size;
	int nnodes = ptree->nlist->size;

	viterbi_parameters_init(&viterbi, nnodes, nrows);

	int i;
	for (i = 0; i < param->num_threads; i++)
	{
		ptparam[i].pfamily = param->pfamily;
		ptparam[i].pcafe = param->pcafe;
		ptparam[i].num_threads = param->num_threads;

		ptparam[i].num_random_samples = globals.num_random_samples;
		ptparam[i].viterbi = &viterbi;
		ptparam[i].pvalue = param->pvalue;
		ptparam[i].from = i;
		ptparam[i].pCD = pCD;
	}

    // Remove threading until process can be more carefully analyzed
    for (i = 0; i < param->num_threads; ++i)
    {
        __cafe_viterbi_thread_func(ptparam + i);
    }

    for (i = 0; i < ptree->nlist->size - 1; i++)
	{
		viterbi.averageExpansion[i] /= param->pfamily->flist->size;
	}

	memory_free(ptparam);
	ptparam = NULL;
}

void cafe_viterbi_print(pCafeParam param, viterbi_parameters& viterbi)
{
	int i, j;
	int size = param->pfamily->flist->size;
	pCafeTree pcafe = param->pcafe;
	pTree ptree = (pTree)pcafe;
	for (i = 0; i < size; i++)
	{
    pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
    cafe_family_set_size(param->pfamily, pitem, pcafe);
		for (j = 1; j < ptree->nlist->size; j += 2)
		{
			pCafeNode pcnode = (pCafeNode)ptree->nlist->array[j];
			pcnode->familysize = viterbi.viterbiNodeFamilysizes[viterbi_parameters::NodeFamilyKey(pcnode->super.super.id, pitem)];
		}
		cafe_tree_string_print(pcafe);
	}
}

/*******************************************************************************
*	Viterbi
*******************************************************************************/
/* this is for finding the ML path
instead of summing the likelihood across all possible innernodes
choose the maximum likelihood across all possible innernodes */
void __cafe_tree_node_compute_viterbi(pTree ptree, pTreeNode ptnode, va_list ap1)
{
  pCafeTree pcafe = (pCafeTree)ptree;
  pCafeNode pcnode = (pCafeNode)ptnode;
  //double lambda = pcafe->lambda;

  double *tree_factors[2];
  tree_factors[0] = (double *)memory_new(pcafe->size_of_factor, sizeof(double));
  tree_factors[1] = (double *)memory_new(pcafe->size_of_factor, sizeof(double));
  int size;
  int s, c, i, j;
  int* rootfamilysizes;
  int* familysizes;

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
    double* factors = tree_factors[0];
    memset(factors, 0, pcafe->size_of_factor * sizeof(double));
    for (s = rootfamilysizes[0], i = 0; s <= rootfamilysizes[1]; s++, i++)
    {
      if (pcnode->familysize < 0) {
        //fprintf(stderr, "family size not set\n");
        pcnode->likelihoods[i] = 1;
        double tmp = 0;
        for (c = pcafe->familysizes[0], j = 0; c <= pcafe->familysizes[1]; c++, j++)
        {
          tmp = square_matrix_get(pcnode->birthdeath_matrix, s, c);
          if (tmp > factors[i])
          {
            factors[i] = tmp;
            pcnode->viterbi[i] = j;
          }
        }
      }
      else {
        if (pcnode->errormodel) {
          memset((void*)pcnode->likelihoods, 0, pcafe->size_of_factor * sizeof(double));
          //int start = pcnode->familysize+pcnode->errormodel->fromdiff;
          //int end = pcnode->familysize+pcnode->errormodel->todiff;
          //if (start<0) {i=0;j=i-start;} else {i=start;j=0;}
          for (j = 0; j<pcafe->size_of_factor; j++) {
            // conditional probability of measuring i=familysize when true count is j
            pcnode->likelihoods[j] = pcnode->errormodel->errormatrix[pcnode->familysize][j];
          }
        }
        else {
          //pcnode->likelihoods[i] = pcnode->bd[s][pcnode->familysize];
          memset((void*)pcnode->likelihoods, 0, pcafe->size_of_factor * sizeof(double));
          pcnode->likelihoods[pcnode->familysize] = 1;
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
    for (idx = 0; idx < 2; idx++)
    {
      {
        factors[idx] = tree_factors[idx];
        memset(factors[idx], 0, pcafe->size_of_factor * sizeof(double));
        for (s = rootfamilysizes[0], i = 0; s <= rootfamilysizes[1]; s++, i++)
        {
          double tmp = 0;
          for (c = familysizes[0], j = 0; c <= familysizes[1]; c++, j++)
          {
            tmp = square_matrix_get(child[idx]->birthdeath_matrix, s, c) * child[idx]->likelihoods[j];
            if (tmp > factors[idx][i])
            {
              factors[idx][i] = tmp;
              child[idx]->viterbi[i] = j;
            }
          }
        }
      }
    }
    size = rootfamilysizes[1] - rootfamilysizes[0] + 1;
    for (i = 0; i < size; i++)
    {
      pcnode->likelihoods[i] = factors[0][i] * factors[1][i];
    }
    //		printf("%s : %d, %d\n", pcnode->super.name, pcnode->super.branchlength, i );
  }

  memory_free(tree_factors[0]);
  memory_free(tree_factors[1]);
}

void __cafe_tree_node_backtrack_viterbi(pTree ptree, pTreeNode ptnode, va_list ap1)
{
  pCafeTree pcafe = (pCafeTree)ptree;
  pCafeNode pcnode = (pCafeNode)ptnode;
  if (tree_is_leaf(ptnode) && (pcnode->familysize >= 0)) return;

  if (ptree->root == ptnode)
  {
    pcnode->familysize = pcafe->rootfamilysizes[0] + __maxidx(pcnode->likelihoods, pcafe->rfsize);
    /* check family size for all nodes first */
    if (pcnode->familysize>10000) {
      fprintf(stderr, "ERROR: FamilySize larger than bd array size Something wrong\n");
      exit(-1);
    }
    /* end check family size for all nodes first */

  }
  else
  {
    pCafeNode pcparent = (pCafeNode)ptnode->parent;
    int base = tree_is_root(ptree, (pTreeNode)pcparent) ? pcafe->rootfamilysizes[0] : pcafe->familysizes[0];
    pcnode->familysize = pcnode->viterbi[pcparent->familysize - base] + pcafe->familysizes[0];
    /* check family size for all nodes first */
    if (pcnode->familysize>10000) {
      fprintf(stderr, "ERROR: FamilySize larger than bd array size Something wrong\n");
      exit(-1);
    }
    /* end check family size for all nodes first */
  }
}

void __cafe_tree_node_compute_clustered_viterbi(pTree ptree, pTreeNode ptnode, va_list ap1)
{
  va_list ap;
  va_copy(ap, ap1);
  int num_likelihoods = va_arg(ap, int);
  va_end(ap);

  pCafeTree pcafe = (pCafeTree)ptree;
  double *tree_factors[2];
  tree_factors[0] = (double *)memory_new(pcafe->size_of_factor, sizeof(double));
  tree_factors[1] = (double *)memory_new(pcafe->size_of_factor, sizeof(double));

  pCafeNode pcnode = (pCafeNode)ptnode;

  int size;
  int i, k;
  int* rootfamilysizes;
  int* familysizes;

  int maxFamilySize = MAX(pcafe->rootfamilysizes[1], pcafe->familysizes[1]);
  if (!chooseln_is_init2(&cache))
  {
    chooseln_cache_init2(&cache, maxFamilySize);
  }
  else if (get_chooseln_cache_size2(&cache) < maxFamilySize)
  {
    chooseln_cache_resize2(&cache, maxFamilySize);
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

    int range = rootfamilysizes[1] - rootfamilysizes[0] + 1;
    initialize_leaf_likelihoods_for_viterbi(pcnode->k_likelihoods, num_likelihoods, range, pcnode->familysize, pcafe->size_of_factor, pcnode->errormodel);
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
        if (k == 0) {
          factors[idx] = tree_factors[idx];
          memset(factors[idx], 0, pcafe->size_of_factor * sizeof(double));
        }
        compute_viterbis(child[idx], k, factors[idx], rootfamilysizes[0], rootfamilysizes[1], familysizes[0], familysizes[1]);
      }
      size = rootfamilysizes[1] - rootfamilysizes[0] + 1;
      for (i = 0; i < size; i++)
      {
        pcnode->k_likelihoods[k][i] = factors[0][i] * factors[1][i];
      }
    }
  }
}


void cafe_tree_viterbi(pCafeTree pcafe)
{
  if (pcafe->super.postfix)
  {
    pArrayList postfix = pcafe->super.postfix;
    pArrayList prefix = pcafe->super.prefix;
    int i = 0;
    for (i = 0; i < postfix->size; i++)
    {
      __cafe_tree_node_compute_viterbi((pTree)pcafe, (pTreeNode)postfix->array[i], NULL);
    }
    for (i = 0; i < prefix->size; i++)
    {
      __cafe_tree_node_backtrack_viterbi((pTree)pcafe, (pTreeNode)prefix->array[i], NULL);
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
  if (pcafe->super.postfix)
  {
    pArrayList postfix = pcafe->super.postfix;
    pArrayList prefix = pcafe->super.prefix;
    int i = 0;
    for (i = 0; i < postfix->size; i++)
    {
      __cafe_tree_node_compute_viterbi((pTree)pcafe, (pTreeNode)postfix->array[i], NULL);
    }
    if (param->posterior) {
      pCafeNode root = (pCafeNode)postfix->array[postfix->size - 1];
      if (tree_is_root((pTree)pcafe, (pTreeNode)root)) {
        for (i = 0; i < param->pcafe->rfsize; i++)	// j: root family size
        {
          // likelihood and posterior both starts from 1 instead of 0 
          root->likelihoods[i] = exp(log(root->likelihoods[i]) + log(param->prior_rfsize[i]));	//prior_rfsize also starts from 1
        }
      }
    }
    for (i = 0; i < prefix->size; i++)
    {
      __cafe_tree_node_backtrack_viterbi((pTree)pcafe, (pTreeNode)prefix->array[i], NULL);
    }
  }
  else
  {
    tree_traveral_postfix((pTree)pcafe, __cafe_tree_node_compute_viterbi);
    tree_traveral_prefix((pTree)pcafe, __cafe_tree_node_backtrack_viterbi);
  }

}

void viterbi_parameters::set_node_familysize(pCafeTree tree, pCafeFamilyItem pItem)
{
  int nnodes = (tree->super.nlist->size - 1) / 2;
  for (int j = 0; j < nnodes; j++)
  {
    pCafeNode pnode = (pCafeNode)tree->super.nlist->array[2 * j + 1];
    pnode->familysize = viterbiNodeFamilysizes[NodeFamilyKey(pnode->super.super.id, pItem)];
  }
}

void viterbi_parameters::compute_size_deltas(pTree ptree, pCafeFamilyItem pitem)
{
  int nnodes = (ptree->nlist->size - 1) / 2;

  for (int j = 0; j < nnodes; j++)
  {
    pCafeNode pcnode = (pCafeNode)ptree->nlist->array[2 * j + 1];
    viterbiNodeFamilysizes[NodeFamilyKey(pcnode->super.super.id, pitem)] = pcnode->familysize;
    pCafeNode child[2] = { (pCafeNode)((pTreeNode)pcnode)->children->head->data,
      (pCafeNode)((pTreeNode)pcnode)->children->tail->data };
    for (int k = 0; k < 2; k++)
    {
      int m = j * 2 + k;
      if (child[k]->familysize > pcnode->familysize) 
        expandRemainDecrease[m].expand++;
      else if (child[k]->familysize == pcnode->familysize) 
        expandRemainDecrease[m].remain++;
      else 
        expandRemainDecrease[m].decrease++;

      pthread_mutex_lock(&mutex_cafe_viterbi);
      averageExpansion[m] += child[k]->familysize - pcnode->familysize;
      pthread_mutex_unlock(&mutex_cafe_viterbi);
    }
  }
}

void viterbi_parameters::clear(int nnodes)
{
  if (!viterbiPvalues.empty())
  {
    averageExpansion.clear();
    if (maximumPvalues)
    {
      memory_free(maximumPvalues);
      maximumPvalues = NULL;
    }
  }
  if (cutPvalues)
  {
    memory_free_2dim((void**)cutPvalues, nnodes, 0, NULL);
  }
  viterbiPvalues.clear();
  cutPvalues = NULL;
  maximumPvalues = NULL;
  expandRemainDecrease.clear();
}

void viterbi_family_print(pCafeTree pcafe, pCafeFamily pfamily, int idx)
{
    cafe_family_set_size_with_family_forced(pfamily, idx, pcafe);
    compute_tree_likelihoods(pcafe);
    int ridx = __maxidx(((pCafeNode)pcafe->super.root)->likelihoods, pcafe->rfsize) + pcafe->rootfamilysizes[0];
    double mlh = __max(((pCafeNode)pcafe->super.root)->likelihoods, pcafe->rfsize);
    cafe_tree_viterbi(pcafe);
    pString pstr = cafe_tree_string(pcafe);
    printf("%g(%d)\t%s\n", mlh, ridx, pstr->buf);
    string_free(pstr);
}

