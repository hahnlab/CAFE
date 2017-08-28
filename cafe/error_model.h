#ifndef ERROR_MODEL_H_D4D61013_69C7_4426_B2D4_84389690EA75
#define ERROR_MODEL_H_D4D61013_69C7_4426_B2D4_84389690EA75

#include <string>
#include <iosfwd>

extern "C" {
#include "cafe.h"
#include "cafe_shell.h"
}

pErrorStruct get_error_model(pCafeFamily family, std::string filename);
int set_error_matrix_from_file(pCafeFamily family, pCafeTree pTree, family_size_range& range, std::string filename, std::string speciesname);
int remove_error_model(pCafeFamily family, pCafeTree pcafe, std::string species_name);
void free_error_model(pCafeFamily family, pCafeTree pcafe);
void init_error_ptr(pCafeFamily family, pCafeTree pTree, pErrorStruct errormodel, std::string speciesname);

std::ostream& operator<<(std::ostream& ost, ErrorStruct& errormodel);
std::istream& operator>>(std::istream& ifst, ErrorStruct& errormodel);

pErrorMeasure estimate_error_double_measure(std::ostream& log, const char* error1, const char* error2, int b_symmetric, int max_diff, int b_peakzero, int maxFamilySize);
pErrorMeasure estimate_error_true_measure(std::ostream& log, const char* errorfile, const char* truefile, int b_symmetric, int max_diff, int b_peakzero, int max_family_size);
int read_freq_from_measures(std::istream* ist1, std::istream* ist2, int* sizeFreq, int maxFamilySize);


#endif
