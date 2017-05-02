#ifndef GENE_FAMILY_H_4DEBB8E1_1D4E_4AA9_8B27_4D8FD59507A6
#define GENE_FAMILY_H_4DEBB8E1_1D4E_4AA9_8B27_4D8FD59507A6

#include <vector>
#include <string>

extern "C" {
#include "cafe.h"
#include "family.h"
}

std::vector<std::string> tokenize(std::string s);

pCafeFamily cafe_family_init(const std::vector<std::string>& species_list);
pCafeFamily load_gene_families(std::istream& ist, int bpatcheck);
void cafe_family_free(pCafeFamily pcf);

void cafe_family_add_item(pCafeFamily pcf, const std::vector<std::string>& data);
void cafe_family_item_free(pCafeFamilyItem pitem);

double cross_validate_by_family(const char* queryfile, const char* truthfile, const char* errortype);

#endif
