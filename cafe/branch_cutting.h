#ifndef BRANCH_CUTTING_H_043A6EEC_78B6_4760_A3B8_F9343BB71F0E
#define BRANCH_CUTTING_H_043A6EEC_78B6_4760_A3B8_F9343BB71F0E

#include <vector>

extern "C" {
#include <family.h>
}

class viterbi_parameters;
class Globals;

typedef std::vector<std::vector<double> > matrix;

class CutBranch
{
	int node_count;
public:
	std::vector<std::pair<matrix, matrix> > pCDSs;

	CutBranch(int nnodes) : node_count(nnodes), pCDSs(nnodes)
	{
	}
};

void set_size_for_split(pCafeFamily pcf, int idx, pCafeTree pcafe);
void cafe_branch_cutting(Globals& globals, viterbi_parameters& viterbi);
void cut_branch(CutBranch& cb, pTree ptree, pCafeTree paramCafe, family_size_range& range, int num_threads, int num_random_samples, int b, std::ostream& ost);
void compute_cutpvalues(pCafeTree pparamcafe, pCafeFamily family, int num_random_samples, int b, int range_start, int range_stop, viterbi_parameters& viterbi, double pvalue, std::vector<double>& p1, double** p2, CutBranch& cb);

#endif
