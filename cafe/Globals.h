#ifndef GLOBALS_H_38BD6CF7_3019_4EA7_9338_FAB9AA0B1B51
#define GLOBALS_H_38BD6CF7_3019_4EA7_9338_FAB9AA0B1B51

extern "C" {
#include "family.h"
}

class viterbi_parameters;
class cross_validator;

class Globals
{
public:
	CafeParam param;

	Globals();
	~Globals();

	pTree mu_tree;

	viterbi_parameters* viterbi;
  cross_validator* validator;

	void Clear(int btree_skip);
	void Prepare();

	int  num_random_samples;

};

#endif
