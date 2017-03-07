#ifndef SIMERROR_H_64CD5F45_CBC7_4A50_868F_6DF76C7E4A31
#define SIMERROR_H_64CD5F45_CBC7_4A50_868F_6DF76C7E4A31

#include <string>
#include <iosfwd>

extern "C" {
#include "cafe.h"
}

void simulate_misclassification(pCafeFamily pcf);
double simerror(pCafeFamily pcf, std::string prefix, int repeat);
size_t get_random(std::vector<double> misclassification_probability);
void write_species_counts(pCafeFamily pcf, std::ostream& ost);

#endif
