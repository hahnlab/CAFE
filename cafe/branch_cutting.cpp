#include <sstream>
#include <utility>

#include "branch_cutting.h"
#include "conditional_distribution.h"

extern "C" {
#include <family.h>
#include "cafe.h"
}

/*
* cdlen: length of conddist : number of trials
*/
double** p_values_of_two_trees(pCafeTree pcafe1, pCafeTree pcafe2,
	double** pvalues, const std::pair<matrix, matrix>& cond_dist,
	int cdlen)
{
	compute_tree_likelihoods(pcafe1);
	compute_tree_likelihoods(pcafe2);
	double* lh1 = get_likelihoods(pcafe1);
	double* lh2 = get_likelihoods(pcafe2);
	int s1, s2, t;
	double p;
	for (s2 = 0; s2 < pcafe1->rfsize; s2++)
	{
		for (s1 = 0; s1 < pcafe1->rfsize; s1++)
		{
			p = 0;
			for (t = 0; t < cdlen; t++)
			{
				p += pvalue(lh1[s1] * lh2[s2] / (cond_dist.second[s2])[t], &cond_dist.first[s1][0], cdlen);
			}
			pvalues[s1][s2] = p / cdlen;
		}
	}
	return pvalues;
}



/**************************************************************************
* BranchCutting
**************************************************************************/
typedef struct
{
	pCafeParam cafeparam;
	CutBranch *cb;
	int range[2];
}BranchCuttingParam;

typedef BranchCuttingParam* pBranchCuttingParam;

void set_size_for_split(pCafeFamily pcf, int idx, pCafeTree pcafe)
{
	int i, j;
	pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[idx];
	for (i = 0; i < pcf->num_species; i++)
	{
		for (j = 0; j < pcafe->super.nlist->size; j += 2)
		{
			pPhylogenyNode pnode = (pPhylogenyNode)pcafe->super.nlist->array[j];
			if (pnode->name[0] & 0x80) continue;
			if (string_pchar_cmp_ignore_case(pnode->name, pcf->species[i]))
			{
				pnode->name[0] |= 0x80;
				((pCafeNode)pnode)->familysize = pitem->count[i];
				break;
			}
		}
	}

	for (j = 0; j < pcafe->super.nlist->size; j += 2)
	{
		pPhylogenyNode pnode = (pPhylogenyNode)pcafe->super.nlist->array[j];
		pnode->name[0] &= 0x7F;
	}
}

pArrayList to_arraylist(matrix& v)
{
	pArrayList result = arraylist_new(10);
	for (size_t i = 0; i < v.size(); ++i)
	{
		double * temp = (double *)calloc(v[i].size(), sizeof(double));
		std::copy(temp, temp + v[i].size(), v[i].begin());
		arraylist_add(result,temp);
	}
	return result;
}

void compute_cutpvalues(pCafeTree pparamcafe, pCafeFamily family, int num_random_samples, int b, int range_start, int range_stop, viterbi_parameters& viterbi, double pvalue, double *p1, double** p2, CutBranch& cb)
{
	pTree ptree = (pTree)pparamcafe;
	if (tree_is_root(ptree, (pTreeNode)ptree->nlist->array[b]))
	{
		return;
	}

	pCafeTree pcafe = cafe_tree_copy(pparamcafe);
	pCafeTree psub = cafe_tree_split(pcafe, b);

	for (int i = range_start; i < range_stop; i++)
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)family->flist->array[i];
		if (pitem->ref >= 0 && pitem->ref != i) continue;
		if (viterbi.maximumPvalues[i] > pvalue)
		{
			viterbi.cutPvalues[b][i] = -1;
			continue;
		}
		if (tree_is_leaf(psub->super.root) || tree_is_leaf(pcafe->super.root))
		{
			pCafeTree pct = tree_is_leaf(psub->super.root) ? pcafe : psub;
			set_size_for_split(family, i, pct);
			pArrayList arr = to_arraylist(cb.pCDSs[b].first);
			assert(cb.pCDSs[b].first.size() == (size_t)pct->rfsize);
			cafe_tree_p_values(pct, p1, arr, num_random_samples);
			arraylist_free(arr, free);
			viterbi.cutPvalues[b][i] = __max(p1, pcafe->rfsize);
		}
		else
		{
			set_size_for_split(family, i, pcafe);
			set_size_for_split(family, i, psub);
			p_values_of_two_trees(pcafe, psub, p2, cb.pCDSs[b], num_random_samples / 10);
			viterbi.cutPvalues[b][i] = 0;
			int m, n;
			for (m = 0; m < pcafe->rfsize; m++)
			{
				for (n = 0; n < pcafe->rfsize; n++)
				{
					if (p2[m][n] > viterbi.cutPvalues[b][i])
					{
						viterbi.cutPvalues[b][i] = p2[m][n];
					}
				}
			}
		}
	}
	cafe_tree_free(pcafe);
	cafe_tree_free(psub);
}

void* __cafe_branch_cutting_thread_func(void* ptr)
{
	pBranchCuttingParam pbc = (pBranchCuttingParam)ptr;
	pCafeParam param = pbc->cafeparam;

#ifdef VERBOSE
	printf("Branch cutting : %d ~ %d\n", pbc->range[0], pbc->range[1] - 1);
#endif

	pTree ptree = (pTree)param->pcafe;
	int nnodes = ptree->nlist->size;
	double* p1 = (double*)memory_new(param->pcafe->rfsize, sizeof(double));
	double** p2 = (double**)memory_new_2dim(param->pcafe->rfsize, param->pcafe->rfsize, sizeof(double));

	for (int b = 0; b < nnodes; b++)
	{
		compute_cutpvalues(param->pcafe, param->pfamily, param->num_random_samples, b, pbc->range[0], pbc->range[1], param->viterbi, param->pvalue, p1, p2, *pbc->cb);
	}
	memory_free(p1);
	p1 = NULL;
	memory_free_2dim((void**)p2, param->pcafe->rfsize, param->pcafe->rfsize, NULL);
	return (NULL);
}

std::ostream& operator<<(std::ostream& os, CafeTree& tree)
{
	pString pstr = cafe_tree_string(&tree);
	os << pstr->buf;
	string_free(pstr);
	return os;
}


void cut_branch(CutBranch& cb, pTree ptree, pCafeTree paramCafe, family_size_range& range, int num_threads, int num_random_samples, int b, std::ostream& ost)
{
	if (tree_is_root(ptree, (pTreeNode)ptree->nlist->array[b]))
	{
		cb.pCDSs[b].first.clear();
		cb.pCDSs[b].second.clear();
		return;
	}
	pCafeTree pcafe = cafe_tree_copy(paramCafe);
	pCafeTree psub = cafe_tree_split(pcafe, b);

	ost << ">> " << b << "  --------------------\n";
	ost << *pcafe << "\n";
	ost << *psub << "\n";

	if (tree_is_leaf(psub->super.root))
	{
		cb.pCDSs[b].first = cafe_conditional_distribution(pcafe, &range, num_threads, num_random_samples);
		cb.pCDSs[b].second.clear();
	}
	else if (tree_is_leaf(pcafe->super.root))
	{
		cb.pCDSs[b].first = cafe_conditional_distribution(psub, &range, num_threads, num_random_samples);
		cb.pCDSs[b].second.clear();
	}
	else
	{
		num_random_samples /= 10;
		cb.pCDSs[b].first = cafe_conditional_distribution(pcafe, &range, num_threads, num_random_samples);
		cb.pCDSs[b].second = cafe_conditional_distribution(psub, &range, num_threads, num_random_samples);
	}

	cafe_tree_free(pcafe);
	cafe_tree_free(psub);
}

void cafe_branch_cutting(pCafeParam param)
{
	cafe_log(param, "Running Branch Cutting....\n");

	pTree ptree = (pTree)param->pcafe;
	int i, b;
	int nnodes = ptree->nlist->size;
	CutBranch cb(nnodes);

	for (b = 0; b < nnodes; b++)
	{
		std::ostringstream ost;
		cut_branch(cb, ptree, param->pcafe, param->family_size, param->num_threads, param->num_random_samples, b, ost);
		cafe_log(param, ost.str().c_str());
	}

	int threadstep = param->pfamily->flist->size / param->num_threads;
	pBranchCuttingParam ptparam = (pBranchCuttingParam)calloc(param->num_threads, sizeof(BranchCuttingParam));

	int nrows = param->pfamily->flist->size;
	param->viterbi.cutPvalues = (double**)memory_new_2dim(nnodes, nrows, sizeof(double));
	int rid = ptree->root->id;
	for (i = 0; i < nrows; i++)
	{
		param->viterbi.cutPvalues[rid][i] = -1;
	}

	int r = 0;
	for (i = 0; i < param->num_threads; i++, r += threadstep)
	{
		ptparam[i].cafeparam = param;
		ptparam[i].cb = &cb;
		ptparam[i].range[0] = r;
		ptparam[i].range[1] = r + threadstep;
	}
	ptparam[param->num_threads - 1].range[1] = param->pfamily->flist->size;
	thread_run(param->num_threads, __cafe_branch_cutting_thread_func, ptparam, sizeof(BranchCuttingParam));

	for (i = 0; i < nrows; i++)
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
		if (pitem->ref < 0 || pitem->ref == i) continue;
		for (b = 0; b < nnodes; b++)
		{
			param->viterbi.cutPvalues[b][i] = param->viterbi.cutPvalues[b][pitem->ref];
		}
	}

	memory_free(ptparam);
	ptparam = NULL;
	cafe_log(param, "Done : Branch Cutting\n");
}

