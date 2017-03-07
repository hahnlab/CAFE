#ifndef CONDITIONAL_DISTRIBUTION_F12ED8EA_52AE_4508_A175_5DDCBE682A46
#define CONDITIONAL_DISTRIBUTION_F12ED8EA_52AE_4508_A175_5DDCBE682A46

#include <vector>

extern "C" {
#include "cafe.h"
#include <family.h>
}

typedef std::vector<std::vector<double> > matrix;

std::vector<double> get_random_probabilities(pCafeTree pcafe, int rootFamilysize, int trials);
matrix conditional_distribution(pCafeTree pcafe, int range_start, int range_end, int num_trials);
matrix cafe_conditional_distribution(pCafeTree pTree, family_size_range *range, int numthreads, int num_random_samples);

#endif
