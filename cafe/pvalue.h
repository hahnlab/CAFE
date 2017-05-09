#ifndef PVALUE_H_B9D3D490_57EB_4E47_8057_0F8D7098FAA4
#define PVALUE_H_B9D3D490_57EB_4E47_8057_0F8D7098FAA4

#include <iosfwd>
#include <vector>

extern "C" {
#include <family.h>
}

void check_cache_and_compute_likelihoods(pCafeTree pTree, int max);
void print_pvalues(std::ostream& ost, pCafeTree pcafe, int max, int num_random_samples);
void read_pvalues(std::istream& ist, int count);
void write_pvalues(std::ostream& ost, pArrayList values, int count);
void pvalues_for_family(pCafeTree pTree, pCafeFamily family, family_size_range *range, int numthreads, int num_random_samples, int index);
void cafe_tree_p_values(pCafeTree pcafe, double* p, pArrayList pconddist, int cdlen);

class ConditionalDistribution
{
public:
	static std::vector<std::vector<double> > matrix;
	static void reset(pCafeTree pTree, family_size_range *range, int numthreads, int num_random_samples);
	static pArrayList to_arraylist();
};

#endif
