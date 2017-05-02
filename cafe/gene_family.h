#ifndef GENE_FAMILY_H_4DEBB8E1_1D4E_4AA9_8B27_4D8FD59507A6
#define GENE_FAMILY_H_4DEBB8E1_1D4E_4AA9_8B27_4D8FD59507A6

extern "C" {
#include "family.h"
}

pCafeFamily cafe_family_init(pArrayList data);
pCafeFamily cafe_family_new(const char* file, int bpatcheck);
double cross_validate_by_family(const char* queryfile, const char* truthfile, const char* errortype);

#endif
