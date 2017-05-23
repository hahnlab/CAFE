#include <vector>

#include "likelihood_ratio.h"

extern "C" {
#include "cafe.h"
#include <family.h>

	extern pBirthDeathCacheArray probability_cache;
}

const double bl_augment = 0.5;


/**************************************************************************
* Likelihood ratio test with more than one lambda
**************************************************************************/

struct LRT2LParam
{
	pCafeParam cafeparam;
	std::vector<int>& lambda;
	std::vector<double> & pvalues;
	param_func lfunc;
	pTree lambda_tree;
	int    num_lambdas;
	std::vector<double*> &lambda_cache;
	pBirthDeathCacheArray* PBDC;
};

void likelihood_ratio_report(pCafeFamily pfamily, 
	pCafeTree pcafe, 
	const std::vector<double> &pvalues, 
	const std::vector<int> &plambda, 
	const std::vector<double*> &lambda_cache, 
	FILE *fout)
{
	int i;
	int fsize = pfamily->flist->size;
	for (i = 0; i < fsize; i++)
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)pfamily->flist->array[i];
		cafe_family_set_size(pfamily, pitem, pcafe);
		pString pstr = cafe_tree_string(pcafe);
		fprintf(fout, "%s\t%s\t", pitem->id, pstr->buf);
		double* l = lambda_cache[plambda[i]];
		fprintf(fout, "(%d, %lf,%lf)\t%g\t%f\n", plambda[i], l[0], l[1], pvalues[i], pvalues[1] == 1 ? 1 : 1 - chi2cdf(pvalues[i], 1));
		string_free(pstr);
	}
}

void update_branchlength(pCafeTree pcafe, pTree lambda_tree, double bl_augment, int *old_branchlength, int* t)
// t[0] is for lambda[1]
// t[1] is for lambda[2]
{
	int i;
	if (lambda_tree)
	{
		pArrayList nlist = pcafe->super.nlist;
		pArrayList lambda_nlist = lambda_tree->nlist;
		for (i = 0; i < nlist->size; i++)
		{
			pPhylogenyNode lambda_pnode = (pPhylogenyNode)lambda_nlist->array[i];
			pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
			old_branchlength[i] = pnode->branchlength;
			if (lambda_pnode->taxaid > 0)
			{
				assert(lambda_pnode->taxaid == 1);
				pnode->branchlength += pnode->branchlength * bl_augment * t[lambda_pnode->taxaid - 1];
			}
		}
	}
}


double __cafe_lhr_get_likelihood_for_diff_lambdas(pCafeParam param, int idx, int t, std::vector<double*> &lambda_cache, std::vector<pBirthDeathCacheArray>& PBDC)
{
	update_branchlength(param->pcafe, param->lambda_tree, bl_augment, param->old_branchlength, &t);
	if (lambda_cache[t] == NULL)
	{
		param->lambda = NULL;
		cafe_best_lambda_by_fminsearch(param, param->num_lambdas, 0);
		lambda_cache[t] = param->lambda;
		reset_birthdeath_cache(param->pcafe, 0, &param->family_size);
		PBDC[t] = probability_cache;
	}
	else
	{
		memcpy(param->lambda, lambda_cache[t], sizeof(double)*param->num_lambdas);
		param->param_set_func(param, param->lambda);
		probability_cache = PBDC[t];
		cafe_tree_set_birthdeath(param->pcafe, probability_cache);
	}
	int i;
  pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[idx];
  cafe_family_set_size(param->pfamily, pitem, param->pcafe);
	compute_tree_likelihoods(param->pcafe);
	double mlh = __max(get_likelihoods(param->pcafe), param->pcafe->rfsize);
	pTree ptree = (pTree)param->pcafe;
	pArrayList nlist = ptree->nlist;
	for (i = 0; i < nlist->size; i++)
	{
		pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
		pnode->branchlength = param->old_branchlength[i];
	}
	return mlh;
}


void* __cafe_lhr_for_diff_lambdas_i(pCafeParam param,
	std::vector<int>& lambda,
	std::vector<double> & pvalues,
	param_func lfunc,
	pTree lambda_tree,
	int    num_lambdas,
	std::vector<double*> &lambda_cache,
	std::vector<pBirthDeathCacheArray> &PBDC
)
{
	int i, j;
	pCafeTree pcafe = cafe_tree_copy(param->pcafe);

	pCafeParam cpy_param = cafe_copy_parameters(param);
	cpy_param->num_lambdas = num_lambdas;
	cpy_param->lambda_tree = lambda_tree;
	cpy_param->param_set_func = lfunc;

	int fsize = param->pfamily->flist->size;
	for (i = 0; i < fsize; i += 1)
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
		if (pitem->ref >= 0 && pitem->ref != i) continue;
		cafe_family_set_size(param->pfamily, pitem, pcafe);
		compute_tree_likelihoods(pcafe);
		double maxlh1 = __max(get_likelihoods(pcafe), pcafe->rfsize);
		double prev = -1;
		double next = __cafe_lhr_get_likelihood_for_diff_lambdas(cpy_param, i, 0, lambda_cache, PBDC);
		for (j = 1; prev < next; j++)
		{
			prev = next;
			next = __cafe_lhr_get_likelihood_for_diff_lambdas(cpy_param, i, j, lambda_cache, PBDC);
		}
		pvalues[i] = (prev == maxlh1) ? 1 : 2 * (log(prev) - log(maxlh1));
		lambda[i] = j - 2;
	}
	cafe_free_copy_parameters(cpy_param);
	cafe_tree_free(pcafe);
	return (NULL);
}


void cafe_lhr_for_diff_lambdas(pCafeParam param, pTree lambda_tree2, int num_lambdas, param_func lfunc)
{
	std::vector<double *> lambda_cache(100);

	cafe_log(param, "Running Likelihood Ratio Test 2....\n");
	int i;

	std::vector<pBirthDeathCacheArray> PBDC(100);
	for (i = 0; i < 100; i++)
	{
		lambda_cache[i] = NULL;
		PBDC[i] = NULL;
	}

	int nrows = param->pfamily->flist->size;
	std::vector<double> pvalues(nrows);
	std::vector<int> plambda(nrows);

	param_func old_func = param->param_set_func;
	param->param_set_func = cafe_lambda_set_default;
	__cafe_lhr_for_diff_lambdas_i(param, plambda, pvalues, lfunc, lambda_tree2, num_lambdas, lambda_cache, PBDC);
	param->param_set_func = old_func;

	int fsize = param->pfamily->flist->size;
	for (i = 0; i < fsize; i++)
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
		if (pitem->ref < 0 || pitem->ref == i) continue;
		pvalues[i] = pvalues[pitem->ref];
		plambda[i] = plambda[pitem->ref];
	}

	likelihood_ratio_report(param->pfamily, param->pcafe, pvalues, plambda, lambda_cache, param->fout);

	for (i = 0; i < 100; i++)
	{
		if (lambda_cache[i])
		{
			memory_free(lambda_cache[i]);
			lambda_cache[i] = NULL;
			birthdeath_cache_array_free(PBDC[i]);
		}
	}
}

