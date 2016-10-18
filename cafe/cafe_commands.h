#ifndef CAFE_COMMANDS_H_B2201016_8693_4379_9225_8FFEC3562AFE
#define CAFE_COMMANDS_H_B2201016_8693_4379_9225_8FFEC3562AFE

#include <vector>
#include <string>

extern "C" {
#include "family.h"
#include "cafe_shell.h"
}

int cafe_cmd_source(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_list(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_echo(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_date(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_exit(pCafeParam param, std::vector<std::string> tokens);
int cafe_cmd_report(pCafeParam param, std::vector<std::string> tokens);

int cafe_shell_dispatch_command(char* cmd);
void list_commands(std::ostream& ost);
void get_report_parameters(report_parameters &params, std::vector<std::string> tokens);

std::vector<std::string> tokenize(std::string s);

#endif

