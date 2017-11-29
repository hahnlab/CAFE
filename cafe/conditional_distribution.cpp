#include <algorithm>
#include <iterator>

#include "conditional_distribution.h"

/**************************************************************************
* Conditional Distribution
**************************************************************************/
/* get likelihood values conditioned on the rootFamilysize for a number of randomly generated families */
std::vector<double> get_random_probabilities(pCafeTree pcafe, int rootFamilysize, int trials)
{
	int old_rfsize = pcafe->rfsize;
	int old_rfsizes[2] = { pcafe->range.root_min, pcafe->range.root_max };
	int old_fsizes[2] = { pcafe->range.min, pcafe->range.max };

	pcafe->rfsize = 1;
	pcafe->range.root_min = rootFamilysize;
	pcafe->range.root_max = rootFamilysize;

	int maxFamilySize = MAX(pcafe->range.root_max, pcafe->range.max);

	std::vector<double> probs(trials);

	for (int i = 0; i < trials; i++)
	{
		int max = cafe_tree_random_familysize(pcafe, rootFamilysize, maxFamilySize);
		if (pcafe->super.nlist)
		{
			pcafe->range.max = MIN(max + MAX(50, max / 5), pcafe->range.max);
		}
		compute_tree_likelihoods(pcafe);
		probs[i] = ((pCafeNode)pcafe->super.root)->likelihoods[0];
	}

	pcafe->rfsize = old_rfsize;
	pcafe->range.root_min = old_rfsizes[0];
	pcafe->range.root_max = old_rfsizes[1];
	pcafe->range.min = old_fsizes[0];
	pcafe->range.max = old_fsizes[1];

	std::sort(probs.begin(), probs.end());

	return probs;
}


/* conditional distribution on a range or root sizes */
std::vector<std::vector<double> > conditional_distribution(pCafeTree pcafe, int range_start, int range_end, int num_trials)
{
	std::vector<std::vector<double> > pal;
	for (int s = range_start; s <= range_end; s++)
	{
		std::vector<double> p = get_random_probabilities(pcafe, s, num_trials);
		pal.push_back(p);
	}
	return pal;
}


/**************************************************************************
* Contidional Distribution
**************************************************************************/

typedef struct
{
	pCafeTree pTree;
	int num_random_samples;
	int range[2];
	std::vector<std::vector<double> > pCD;
}CDParam;
typedef CDParam* pCDParam;


void* __cafe_conditional_distribution_thread_func(void* ptr)
{
	pCDParam param = (pCDParam)ptr;
	pCafeTree pcafe = cafe_tree_copy(param->pTree);
#ifdef VERBOSE
	printf("CD: %d ~ %d\n", param->range[0], param->range[1]);
#endif
	param->pCD = conditional_distribution(pcafe, param->range[0], param->range[1], param->num_random_samples);
	cafe_tree_free(pcafe);
	return (NULL);
}

matrix cafe_conditional_distribution(pCafeTree pTree, family_size_range *range, int numthreads, int num_random_samples)
{
	int threadstep = pTree->rfsize / numthreads;
	if (threadstep == 0)
	{
		numthreads = pTree->rfsize;
	}
	else
	{
		threadstep--;
	}

	pCDParam ptparam = (pCDParam)memory_new(numthreads, sizeof(CDParam));
	int i, r = range->root_min;
	for (i = 0; i < numthreads; i++, r += threadstep + 1)
	{
		ptparam[i].pTree = pTree;
		ptparam[i].num_random_samples = num_random_samples;
		ptparam[i].range[0] = r;
		ptparam[i].range[1] = r + threadstep;
	}
	ptparam[numthreads - 1].range[1] = range->root_max;
	thread_run(numthreads, __cafe_conditional_distribution_thread_func, ptparam, sizeof(CDParam));
	matrix cdlist;
	for (i = 0; i < numthreads; i++)
	{
		for (size_t j = 0; j < ptparam[i].pCD.size(); j++)
		{
			cdlist.push_back(ptparam[i].pCD[j]);
		}
	}
	memory_free(ptparam);
	ptparam = NULL;
	return cdlist;
}

