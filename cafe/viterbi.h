#ifndef VITERBI_H_A08989A1_B4B4_461C_B863_A1AE2FE9BD98
#define VITERBI_H_A08989A1_B4B4_461C_B863_A1AE2FE9BD98

#include <vector>

extern "C"
{
#include <family.h>
void compute_viterbis(pCafeNode node, int k, double *factors, int rootfamilysize_start, int rootfamilysize_end, int familysize_start, int familysize_end);
}

struct change
{
	change(int e, int r, int d) : expand(e), remain(r), decrease(d)
	{

	}
	change() : expand(0), remain(0), decrease(0)
	{

	}
	int expand;
	int remain;
	int decrease;
};

class viterbi_parameters
{
public:
	/** Number of nodes in the tree */
	int num_nodes;

	/** Number of gene families for which to keep data */
	int num_rows;

	/** Matrix of calculated P values for each node in the tree and each gene family  */
	double** viterbiPvalues;

	/** array of three integers expand, remain, and decrease for each node in the tree relative to its parent */
	std::vector<change> expandRemainDecrease;

	int** viterbiNodeFamilysizes;
	double* maximumPvalues;
	std::vector<double> averageExpansion;
	double** cutPvalues;
} ;

void viterbi_parameters_init(viterbi_parameters *viterbi, int nnodes, int nrows);

void viterbi_set_max_pvalue(viterbi_parameters* viterbi, int index, double val);
void viterbi_parameters_clear(viterbi_parameters* viterbi, int nnodes);
pArrayList cafe_viterbi(pCafeParam param, viterbi_parameters& viterbi, pArrayList pCD);


#endif
