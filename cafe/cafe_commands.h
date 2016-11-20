#ifndef CAFE_COMMANDS_H_B2201016_8693_4379_9225_8FFEC3562AFE
#define CAFE_COMMANDS_H_B2201016_8693_4379_9225_8FFEC3562AFE

#include <vector>
#include <string>
#include <iosfwd>

extern "C" {
#include "family.h"
}

typedef int(*cafe_command2)(pCafeParam cafe_param, std::vector<std::string>);


int cafe_cmd_source(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_list(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_echo(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_date(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_exit(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_gainloss(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_generate_random_family(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_log(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_version(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_print_param(pCafeParam param, std::vector<std::string> tokens);


int cafe_shell_dispatch_command(char* cmd);
void list_commands(std::ostream& ost);

std::vector<std::string> tokenize(std::string s);

// these functions should be moved to another file
void clear_tree_viterbis(pCafeTree psum);
int get_num_trials(std::vector<std::string> args);
int write_family_gainloss(std::ostream& ofst, std::string family_id, pCafeTree tree1, pCafeTree tree2);
void set_node_familysize(pCafeTree tree, int** node_family_sizes, int i);
std::vector<int> get_clusters(int parameterized_k_value, int num_families, double* k_weights);
void write_node_headers(std::ostream& s1, std::ostream& s2, pCafeTree pcafe);
void write_leaves(std::ostream& ofst, pCafeTree pcafe, int *k, int i, int id, bool evens);
void write_version(std::ostream &ost);

#endif

