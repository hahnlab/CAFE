#ifndef CROSS_VALIDATOR_H_F3DD8C50_D7E9_479F_A19D_4D1CBF4F5190
#define CROSS_VALIDATOR_H_F3DD8C50_D7E9_479F_A19D_4D1CBF4F5190

#include <string>

extern "C" {
#include "cafe.h"
}

class cross_validator
{
  char* cv_species_name;
  pArrayList cv_test_species_list;
  pArrayList cv_test_count_list;

  void read_validate_species(const char* file);
  void read_query_family(const char* file);
public:
  void clean_by_family(std::string file, int cv_fold);
  void clean_by_species(std::string file);
  double validate_by_species(pCafeParam param, const char* validatefile, const char* errortype);
  double validate_by_family(pCafeParam param, const char* queryfile, const char* truthfile, const char* errortype);

  std::string get_species_name() const { return cv_species_name; }
};

#endif
