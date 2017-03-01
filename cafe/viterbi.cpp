#include "viterbi.h"

extern "C" {
#include <family.h>
#include "cafe.h"
#include <pthread.h>
}


/**************************************************************************
* Viterbi
**************************************************************************/
void viterbi_parameters_init(viterbi_parameters *viterbi, int nnodes, int nrows)
{
	viterbi->num_nodes = nnodes;
	viterbi->num_rows = nrows;
	viterbi->viterbiPvalues = (double**)memory_new_2dim(nnodes, nrows, sizeof(double));
	viterbi->viterbiNodeFamilysizes = (int**)memory_new_2dim(nnodes, nrows, sizeof(int));
	viterbi->maximumPvalues = (double*)memory_new(nrows, sizeof(double));
	viterbi->averageExpansion.clear();
	viterbi->expandRemainDecrease.clear();
	viterbi->averageExpansion.resize(nnodes);
	viterbi->expandRemainDecrease.resize(nnodes);
}

void viterbi_parameters_clear(viterbi_parameters* viterbi, int nnodes)
{
	//	viterbi_parameters* viterbi = &param->viterbi;
	if (viterbi->viterbiPvalues)
	{
		int num = (nnodes - 1) / 2;
		memory_free_2dim((void**)viterbi->viterbiPvalues, num, 0, NULL);
		memory_free_2dim((void**)viterbi->viterbiNodeFamilysizes, num, 0, NULL);
		viterbi->averageExpansion.clear();
		if (viterbi->maximumPvalues)
		{
			memory_free(viterbi->maximumPvalues);
			viterbi->maximumPvalues = NULL;
		}
	}
	if (viterbi->cutPvalues)
	{
		memory_free_2dim((void**)viterbi->cutPvalues, nnodes, 0, NULL);
	}
	viterbi->viterbiPvalues = NULL;
	viterbi->viterbiNodeFamilysizes = NULL;
	viterbi->cutPvalues = NULL;
	viterbi->maximumPvalues = NULL;
	viterbi->expandRemainDecrease.clear();
}

void viterbi_set_max_pvalue(viterbi_parameters* viterbi, int index, double val)
{
	assert(index < viterbi->num_rows);
	viterbi->maximumPvalues[index] = val;
}


typedef struct
{
	pCafeFamily pfamily;
	pCafeTree pcafe;
	int num_threads;
	viterbi_parameters *viterbi;
	int num_random_samples;
	double pvalue;

	pArrayList pCD;
	int from;
}ViterbiParam;

typedef ViterbiParam*  pViterbiParam;

pthread_mutex_t mutex_cafe_viterbi = PTHREAD_MUTEX_INITIALIZER;

void viterbi_section(pCafeFamily pcf, double pvalue, int num_random_samples, viterbi_parameters *viterbi, int i, pCafeTree pcafe, double *cP, pArrayList pCD)
{
	pTree ptree = (pTree)pcafe;
	int nnodes = (ptree->nlist->size - 1) / 2;

	cafe_family_set_size_with_family_forced(pcf, i, pcafe);

	cafe_tree_p_values(pcafe, cP, pCD, num_random_samples);
	viterbi_set_max_pvalue(viterbi, i, __max(cP, pcafe->rfsize));
	cafe_tree_viterbi(pcafe);
	/* check family size for all nodes first */
	for (int j = 0; j < nnodes; j++)
	{
		pCafeNode pcnode = (pCafeNode)ptree->nlist->array[j];
		if (pcnode->familysize>10000) {
			fprintf(stderr, "ERROR: FamilySize larger than bd array size Something wrong\n");
			exit(-1);
		}
	}
	/* end check family size for all nodes first */

	for (int j = 0; j < nnodes; j++)
	{
		pCafeNode pcnode = (pCafeNode)ptree->nlist->array[2 * j + 1];
		viterbi->viterbiNodeFamilysizes[j][i] = pcnode->familysize;
		pCafeNode child[2] = { (pCafeNode)((pTreeNode)pcnode)->children->head->data,
			(pCafeNode)((pTreeNode)pcnode)->children->tail->data };
		for (int k = 0; k < 2; k++)
		{
			int m = j * 2 + k;
			if (child[k]->familysize > pcnode->familysize) viterbi->expandRemainDecrease[m].expand++;
			else if (child[k]->familysize == pcnode->familysize) viterbi->expandRemainDecrease[m].remain++;
			else viterbi->expandRemainDecrease[m].decrease++;

			pthread_mutex_lock(&mutex_cafe_viterbi);
			viterbi->averageExpansion[m] += child[k]->familysize - pcnode->familysize;
			pthread_mutex_unlock(&mutex_cafe_viterbi);
		}
	}

	if (viterbi->maximumPvalues[i] > pvalue)
	{
		for (int j = 0; j < ptree->nlist->size - 1; j++)
		{
			viterbi->viterbiPvalues[j][i] = -1;
		}
		return;
	}

	for (int j = 0; j < nnodes; j++)
	{
		pCafeNode pcnode = (pCafeNode)ptree->nlist->array[2 * j + 1];
		pCafeNode child[2] = { (pCafeNode)((pTreeNode)pcnode)->children->head->data,
			(pCafeNode)((pTreeNode)pcnode)->children->tail->data };
		for (int k = 0; k < 2; k++)
		{
			double p = square_matrix_get(child[k]->birthdeath_matrix, pcnode->familysize, child[k]->familysize);
			int n = 2 * j + k;
			for (int m = 0; m <= pcafe->familysizes[1]; m++)
			{
				if (square_matrix_get(child[k]->birthdeath_matrix, pcnode->familysize, m) == p)
				{
					viterbi->viterbiPvalues[n][i] += square_matrix_get(child[k]->birthdeath_matrix, pcnode->familysize, m) / 2.0;
				}
				else if (square_matrix_get(child[k]->birthdeath_matrix, pcnode->familysize, m) < p)
				{
					viterbi->viterbiPvalues[n][i] += square_matrix_get(child[k]->birthdeath_matrix, pcnode->familysize, m);
				}
			}
		}
	}
}


void* __cafe_viterbi_thread_func(void* ptr)
{
	pViterbiParam pv = (pViterbiParam)ptr;

	pArrayList pCD = pv->pCD;
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

pArrayList cafe_viterbi(pCafeParam param, viterbi_parameters& viterbi, pArrayList pCD)
{
	cafe_log(param, "Running Viterbi algorithm....\n");

	pViterbiParam ptparam = (pViterbiParam)memory_new(param->num_threads, sizeof(ViterbiParam));
	pTree ptree = (pTree)param->pcafe;

	int nrows = param->pfamily->flist->size;
	int nnodes = ptree->nlist->size - 1;

	viterbi_parameters_init(&viterbi, nnodes, nrows);

	int i;
	for (i = 0; i < param->num_threads; i++)
	{
		ptparam[i].pfamily = param->pfamily;
		ptparam[i].pcafe = param->pcafe;
		ptparam[i].num_threads = param->num_threads;

		ptparam[i].num_random_samples = param->num_random_samples;
		ptparam[i].viterbi = &viterbi;
		ptparam[i].pvalue = param->pvalue;
		ptparam[i].from = i;
		ptparam[i].pCD = pCD;
	}
	thread_run(param->num_threads, __cafe_viterbi_thread_func, ptparam, sizeof(ViterbiParam));

	for (i = 0; i < ptree->nlist->size - 1; i++)
	{
		viterbi.averageExpansion[i] /= param->pfamily->flist->size;
	}

	memory_free(ptparam);
	ptparam = NULL;
	return pCD;
}

void cafe_viterbi_print(pCafeParam param, viterbi_parameters& viterbi)
{
	int i, j;
	int size = param->pfamily->flist->size;
	pCafeTree pcafe = param->pcafe;
	pTree ptree = (pTree)pcafe;
	for (i = 0; i < size; i++)
	{
		cafe_family_set_size(param->pfamily, i, pcafe);
		for (j = 1; j < ptree->nlist->size; j += 2)
		{
			pCafeNode pcnode = (pCafeNode)ptree->nlist->array[j];
			pcnode->familysize = viterbi.viterbiNodeFamilysizes[j / 2][i];
		}
		cafe_tree_string_print(pcafe);
	}
}
