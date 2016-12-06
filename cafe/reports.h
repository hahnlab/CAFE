#ifndef REPORTS_C6B4B223_5C6D_459E_B2D8_CAC46A40A5B4
#define REPORTS_C6B4B223_5C6D_459E_B2D8_CAC46A40A5B4

#include <vector>
#include <string>

extern "C" {
#include "cafe_shell.h"
#include "family.h"
}

typedef struct
{
	int bc, lh, lh2, just_save;
	std::string name;
} report_parameters;


void get_report_parameters(report_parameters &params, std::vector<std::string> tokens);
int cafe_cmd_report(pCafeParam param, std::vector<std::string> tokens);
void cafe_report(pCafeParam param, std::ostream& report_file);
void write_viterbi(std::ostream& ost, const viterbi_parameters& viterbi);
void write_families_header(std::ostream& ost, double **cutPvalues, double**likelihoodRatios);
void write_families_line(std::ostream& ost, pCafeParam param, int i, std::string node_id);
void cafe_do_report(pCafeParam param, report_parameters* params);


#endif
