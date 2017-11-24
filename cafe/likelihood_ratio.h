#ifndef LIKELIHOOD_RATIO_H_7F6D20CB_1791_422E_B34B_54DF43D11DA5
#define LIKELIHOOD_RATIO_H_7F6D20CB_1791_422E_B34B_54DF43D11DA5

#include <vector>

extern "C" {
#include <family.h>
}
void cafe_lhr_for_diff_lambdas(pCafeParam param, pTree lambda_tree2, int num_lambdas, OPTIMIZER_INIT_TYPE lfunc);
void update_branchlength(pCafeTree pcafe, pTree lambda_tree, double bl_augment, int *old_branchlength, int* t);

void likelihood_ratio_report(pCafeFamily pfamily,
	pCafeTree pcafe,
	const std::vector<double> &pvalues,
	const std::vector<int> &plambda,
	const std::vector<double*> &lambda_cache,
	FILE *fout);



#endif