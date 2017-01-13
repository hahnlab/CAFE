#include "pvalue.h"
#include <istream>
#include <ostream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>

extern "C" {
#include "cafe.h"

	extern pBirthDeathCacheArray probability_cache;
}

pArrayList ConditionalDistribution::cafe_pCD;

void check_cache_and_compute_likelihoods(pCafeTree pTree, int max)
{
	family_size_range range;
	range.min = 0;
	range.root_min = 1;
	range.root_max = max + MAX(50, max / 4);
	range.max = max + MAX(50, max / 5);

	if (probability_cache == NULL || probability_cache->maxFamilysize <  MAX(range.root_max, range.max))
	{
		cafe_tree_set_parameters(pTree, &range, 0);
		if (probability_cache)
		{
			int remaxFamilysize = MAX(range.max, range.root_max);
			birthdeath_cache_resize(probability_cache, remaxFamilysize);
			cafe_tree_set_birthdeath(pTree);
		}
		else
		{
			reset_birthdeath_cache(pTree, 0, &range);
		}
	}
	else
	{
		copy_range_to_tree(pTree, &range);
	}
	compute_tree_likelihoods(pTree);
}


void print_pvalues(std::ostream& ost, pCafeTree pcafe, int max, int num_random_samples)
{
	check_cache_and_compute_likelihoods(pcafe, max);
	double* lh = get_likelihoods(pcafe);
	int rfsize = __maxidx(lh, pcafe->rfsize);
	double mlh = lh[rfsize];
	rfsize += pcafe->rootfamilysizes[0];

	double* pcd = cafe_tree_random_probabilities(pcafe, rfsize, num_random_samples);
	double pv = pvalue(mlh, pcd, num_random_samples);
	pString pstr = cafe_tree_string(pcafe);
	ost << pstr->buf << "\n";
	ost << "Root size: " << rfsize << " with maximum likelihood : " << mlh << "\n";
	ost << "p-value: " << pv << "\n";
	memory_free(pcd);
	pcd = NULL;
	string_free(pstr);
}

void write_pvalues(std::ostream& ost, pArrayList values, int count)
{
	ost << std::setw(10) << std::setprecision(9);
	for (int i = 0; i < values->size; i++)
	{
		double* data = (double*)values->array[i];
		ost << data[0];
		for (int j = 1; j < count; j++)
		{
			ost << "\t" << data[j];
		}
		ost << "\n";
	}

}


void read_pvalues(std::istream& ist, int count)
{
	if (ConditionalDistribution::cafe_pCD) arraylist_free(ConditionalDistribution::cafe_pCD, free);
	ConditionalDistribution::cafe_pCD = arraylist_new(11000);

	std::string str;
	while (std::getline(ist, str))
	{
		double * data = (double*)memory_new(count, sizeof(double));
		std::istringstream iss(str);
		for (int i = 0; i < count; ++i)
			iss >> data[i];
		arraylist_add(ConditionalDistribution::cafe_pCD, data);
	}
}

void pvalues_for_family(pCafeTree pTree, pCafeFamily family, family_size_range *range, int numthreads, int num_random_samples, int index)
{
	if (ConditionalDistribution::cafe_pCD == NULL)
	{
		ConditionalDistribution::cafe_pCD = cafe_conditional_distribution(pTree, range, numthreads, num_random_samples);
	}
	cafe_tree_set_parameters(pTree, range, 0);

	cafe_family_set_size(family, index, pTree);

	compute_tree_likelihoods(pTree);
	double* lh = get_likelihoods(pTree);

	std::vector<double> pvalues(pTree->rfsize);

	for (int s = 0; s < pTree->rfsize; s++)
	{
		pvalues[s] = pvalue(lh[s], (double*)ConditionalDistribution::cafe_pCD->array[s], num_random_samples);
	}

	for (int i = 0; i < pTree->rfsize; i++)
	{
		printf("%d\t%lg\t%lg\n", i + pTree->rootfamilysizes[0], lh[i], pvalues[i]);
	}
}