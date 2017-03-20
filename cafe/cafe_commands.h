#ifndef CAFE_COMMANDS_H_B2201016_8693_4379_9225_8FFEC3562AFE
#define CAFE_COMMANDS_H_B2201016_8693_4379_9225_8FFEC3562AFE

#include <vector>
#include <string>
#include <map>
#include <iosfwd>
#include <stdexcept>

extern "C" {
#include "family.h"
#include "cafe_shell.h"
}

class Globals;

typedef int(*cafe_command2)(Globals& globals, std::vector<std::string>);

#define MAKE_FN_NAME(x) int cafe_cmd_##x (Globals& globals, std::vector<std::string>)
#define COMMAND(signal) MAKE_FN_NAME(signal)

COMMAND(accuracy);
COMMAND(branchlength);
COMMAND(cvfamily);
COMMAND(cvspecies);
COMMAND(date);
COMMAND(echo);
COMMAND(errormodel);
COMMAND(esterror);
COMMAND(exit);
COMMAND(extinct);
COMMAND(family);
COMMAND(gainloss);
COMMAND(generate_random_family);
COMMAND(lambdamu);
COMMAND(lhtest);
COMMAND(list);
COMMAND(load);
COMMAND(log);
COMMAND(noerrormodel);
COMMAND(print_param);
COMMAND(pvalue);
COMMAND(retrieve);
COMMAND(rootdist);
COMMAND(save);
COMMAND(score);
COMMAND(simerror);
COMMAND(simextinct);
COMMAND(source);
COMMAND(tree);
COMMAND(version);
COMMAND(viterbi);

std::map<std::string, cafe_command2> get_dispatcher();
int cafe_shell_dispatch_command(Globals& globals, char* cmd);
void list_commands(std::ostream& ost);

std::vector<std::string> tokenize(std::string s);
std::vector<Argument> build_argument_list(std::vector<std::string> tokens);

// these functions should be moved to another file
void clear_tree_viterbis(pCafeTree psum);
int get_num_trials(std::vector<std::string> args);
int write_family_gainloss(std::ostream& ofst, std::string family_id, pCafeTree tree1, pCafeTree tree2);
void set_node_familysize(pCafeTree tree, int** node_family_sizes, int i);
std::vector<int> get_clusters(int parameterized_k_value, int num_families, double* k_weights);
void write_node_headers(std::ostream& s1, std::ostream& s2, pCafeTree pcafe);
void write_leaves(std::ostream& ofst, pCafeTree pcafe, int *k, int i, int id, bool evens);
void write_version(std::ostream &ost);
void write_family(std::ostream& ost, pCafeFamily family);
void viterbi_write(std::ostream& ost, pCafeTree pcafe, pCafeFamily pfamily);
void tree_set_branch_lengths(pCafeTree pcafe, std::vector<int> lengths);

struct load_args {
	int num_threads;
	int num_random_samples;
	double pvalue;
	bool filter;
	std::string log_file_name;
	std::string family_file_name;
};

struct viterbi_args {
	bool all;
	std::string file;
	int idx;
	std::string item_id;
};

struct pvalue_args
{
	std::string infile;
	std::string outfile;
	int index;
};

struct lhtest_args
{
	std::string directory;
	std::string tree;
	std::string outfile;
	double lambda;
};

struct esterror_args
{
	std::string outfile;
	std::vector<std::string> data_error_files;
	bool symmetric;
	bool peakzero;
	std::string truth_file;
	int max_diff;
};

load_args get_load_arguments(std::vector<Argument> pargs);
viterbi_args get_viterbi_arguments(std::vector<Argument> pargs);
pvalue_args get_pvalue_arguments(std::vector<Argument> pargs);
lhtest_args get_lhtest_arguments(std::vector<Argument> pargs);
esterror_args get_esterror_arguments(std::vector<Argument> pargs);
void validate(esterror_args args);

struct roots
{
	std::vector<double> extinct;
	std::vector<int> size;
	std::vector<int> num;
	std::vector<double> avg_extinct;
	int total_extinct;

	pHistogram* phist_data;
	pHistogram* phist_sim;
};

void run_viterbi_sim(pCafeTree pcafe, pCafeFamily pfamily, roots& roots);
int init_histograms(int rfsize, roots& roots, int nsamples);
void get_doubles_array(std::vector<double>& loc, pArgument parg);

const int REQUIRES_FAMILY = 0x01;
const int REQUIRES_TREE = 0x02;
const int REQUIRES_LAMBDA = 0x04;
const int REQUIRES_ERRORMODEL = 0x08;

void prereqs(pCafeParam param, int flags);

class io_error : public std::runtime_error
{
	std::string text;
public:
	io_error(std::string source, std::string file, bool write);
	virtual ~io_error() throw() {}

	virtual const char* what() const throw()
	{
		return text.c_str();
	}
};


#endif

