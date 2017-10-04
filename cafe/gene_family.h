#ifndef GENE_FAMILY_H_4DEBB8E1_1D4E_4AA9_8B27_4D8FD59507A6
#define GENE_FAMILY_H_4DEBB8E1_1D4E_4AA9_8B27_4D8FD59507A6

#include <vector>
#include <string>
#include <iosfwd>

extern "C" {
#include "cafe.h"
#include "family.h"
}

const int REGULAR_WHITESPACE = 0;
const int COMMA_AS_WHITESPACE = 1;

std::vector<std::string> tokenize(std::string s, int flags);

pCafeFamily cafe_family_init(const std::vector<std::string>& species_list);
pCafeFamily load_gene_families(std::istream& ist, char separator, int max_size);
void cafe_family_free(pCafeFamily pcf);

struct gene_family
{
	std::string id;
	std::string desc;
	std::vector<int> values;

	gene_family()
	{
			
	}

	gene_family(std::string i, std::string d, std::vector<int> v) : id(i), desc(d), values(v)
	{
	}
};

std::istream& operator>>(std::istream& ist, gene_family& fam);

void cafe_family_add_item(pCafeFamily pcf, const gene_family& gf);
void cafe_family_item_free(pCafeFamilyItem pitem);

int log_cluster_membership(pCafeFamily pcf, int k_value, double **p_z_membership, std::ostream& log);

void cafe_family_reset_maxlh(pCafeFamily pcf);
int cafe_family_get_index(pCafeFamily pcf, const char* szid);
#endif
