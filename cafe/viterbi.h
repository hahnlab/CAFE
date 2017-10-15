#ifndef VITERBI_H_A08989A1_B4B4_461C_B863_A1AE2FE9BD98
#define VITERBI_H_A08989A1_B4B4_461C_B863_A1AE2FE9BD98

#include <vector>
#include <map>

extern "C"
{
#include <family.h>
#include <cafe.h>
void compute_viterbis(pCafeNode node, int k, double *factors, int rootfamilysize_start, int rootfamilysize_end, int familysize_start, int familysize_end);
}

class Globals;

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

	/** array of three integers expand, remain, and decrease for each node in the tree relative to its parent */
	std::vector<change> expandRemainDecrease;

	// first value is node ID
  using NodeFamilyKey = std::pair<int, pCafeFamilyItem>;

  std::map<NodeFamilyKey, int> viterbiNodeFamilysizes;

  /** Matrix of calculated P values for each node in the tree and each gene family  */
  std::map<NodeFamilyKey, double> viterbiPvalues;

	double* maximumPvalues;
	std::vector<double> averageExpansion;
	double** cutPvalues;

  void set_node_familysize(pCafeTree tree, pCafeFamilyItem pItem);
  void compute_size_deltas(pTree ptree, pCafeFamilyItem pitem);
  void clear(int nnodes);
} ;

void viterbi_parameters_init(viterbi_parameters *viterbi, int nnodes, int nrows);

void viterbi_set_max_pvalue(viterbi_parameters* viterbi, int index, double val);
pArrayList cafe_viterbi(Globals& globals, viterbi_parameters& viterbi, pArrayList pCD);
void viterbi_sum_probabilities(viterbi_parameters *viterbi, pCafeNode pcnode, pCafeFamilyItem item, int max_family_size);
void* __cafe_viterbi_thread_func(void* ptr);
void viterbi_family_print(pCafeTree pcafe, pCafeFamily pfamily, int idx);

struct ViterbiParam
{
	pCafeFamily pfamily;
	pCafeTree pcafe;
	int num_threads;
	viterbi_parameters *viterbi;
	int num_random_samples;
	double pvalue;

	pArrayList pCD;
	int from;
};


#endif
