#ifndef REPORTS_C6B4B223_5C6D_459E_B2D8_CAC46A40A5B4
#define REPORTS_C6B4B223_5C6D_459E_B2D8_CAC46A40A5B4

#include <vector>
#include <string>

extern "C" {
#include "cafe_shell.h"
#include "family.h"
}

#include "viterbi.h"


class Globals;
class viterbi_parameters;

struct family_line_item
{
	std::string node_id;
	std::string tree;
	double max_p_value;
	std::vector<std::pair<double, double> > pvalues;
	std::vector<double> cut_pvalues;
	std::vector<double> likelihood_ratios;

	family_line_item(pCafeFamily family, pCafeTree pcafe, double** likelihoodRatios, viterbi_parameters& viterbi, int i, std::string node_id);
	family_line_item() : max_p_value(0.0) {}
};

struct Report
{
  enum Formats { Unknown, Text, HTML, JSON };
  static int report_format;

	std::string tree;
	std::string lambda_tree;
	std::string id_tree;
	std::vector<double> lambdas;
	std::vector<double> averageExpansion;
	std::vector<change> changes;
	std::vector<std::pair<int, int> > node_pairs;
	std::vector<family_line_item> family_line_items;
	std::vector<int> branch_cutting_output_format;

	Report() {}
	Report(pCafeParam param, viterbi_parameters& viterbi);
};

struct report_parameters
{
  bool branchcutting;
  bool likelihood;
  bool lh2;
  bool just_save;
  Report::Formats format;
  std::string name;
};


report_parameters get_report_parameters(std::vector<std::string> tokens);
int cafe_cmd_report(Globals& globals, std::vector<std::string> tokens);
void write_viterbi(std::ostream& ost, const Report& viterbi);
void write_families_header(std::ostream& ost, bool cutPvalues, bool likelihoodRatios);
void cafe_do_report(Globals& globals, viterbi_parameters& viterbi, report_parameters* params);
int cafe_report_retrieve_data(const char* file, pCafeParam param, viterbi_parameters& viterbi);

std::ostream& operator<<(std::ostream& ost, const family_line_item& item);
std::ostream& operator<<(std::ostream& ost, const Report& report);

/// These I/O manipulators allow selecting the report format
std::ios_base& json(std::ios_base& os);
std::ios_base& html(std::ios_base& os);
#endif
