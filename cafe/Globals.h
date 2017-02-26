#ifndef GLOBALS_H_38BD6CF7_3019_4EA7_9338_FAB9AA0B1B51
#define GLOBALS_H_38BD6CF7_3019_4EA7_9338_FAB9AA0B1B51

extern "C" {
#include "family.h"
}

class Globals
{
public:
	CafeParam param;

	Globals();

	pTree mu_tree;

	void Clear(int btree_skip);
	void Prepare();

};

#endif
