#ifndef GLOBALS_H_38BD6CF7_3019_4EA7_9338_FAB9AA0B1B51
#define GLOBALS_H_38BD6CF7_3019_4EA7_9338_FAB9AA0B1B51

extern "C" {
#include "family.h"
}

class viterbi_parameters;

class Globals
{
public:
	CafeParam param;

	Globals();
	~Globals();

	pTree mu_tree;

	viterbi_parameters* viterbi;

	void Clear(int btree_skip);
	void Prepare();

	int  num_random_samples;

};

#endif
