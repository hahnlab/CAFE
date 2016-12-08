#ifndef CAFE_COMMANDS_H_B2201016_8693_4379_9225_8FFEC3562AFE
#define CAFE_COMMANDS_H_B2201016_8693_4379_9225_8FFEC3562AFE

#include <vector>
#include <string>
#include <iosfwd>

extern "C" {
#include "family.h"
#include "cafe_shell.h"
}

typedef int(*cafe_command2)(pCafeParam cafe_param, std::vector<std::string>);

#define MAKE_FN_NAME(x) int cafe_cmd_##x (pCafeParam cafe_param, std::vector<std::string>)
#define COMMAND(signal) MAKE_FN_NAME(signal)

COMMAND(source);
COMMAND(list);
COMMAND(date);
COMMAND(echo);
COMMAND(exit);
COMMAND(gainloss);
COMMAND(generate_random_family);
COMMAND(log);
COMMAND(version);
COMMAND(print_param);
COMMAND(load);
COMMAND(family);
COMMAND(save);
COMMAND(tree);
COMMAND(score);
COMMAND(viterbi);

int cafe_shell_dispatch_command(char* cmd);
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

struct load_args get_load_arguments(std::vector<Argument> pargs);
struct viterbi_args get_viterbi_arguments(std::vector<Argument> pargs);

#endif

