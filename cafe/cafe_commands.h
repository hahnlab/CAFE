#ifndef CAFE_COMMANDS_H_B2201016_8693_4379_9225_8FFEC3562AFE
#define CAFE_COMMANDS_H_B2201016_8693_4379_9225_8FFEC3562AFE

#include <vector>
#include <string>
#include <iosfwd>

extern "C" {
#include "family.h"
}

int cafe_cmd_source(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_list(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_echo(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_date(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_exit(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_gainloss(pCafeParam param, std::vector<std::string> tokens);

int cafe_shell_dispatch_command(char* cmd);
void list_commands(std::ostream& ost);

std::vector<std::string> tokenize(std::string s);

void clear_tree_viterbis(pCafeTree psum);

int write_family_gainloss(std::ostream& ofst, std::string family_id, pCafeTree tree1, pCafeTree tree2);
void set_node_familysize(pCafeTree tree, int** node_family_sizes, int i);

#endif

