#ifndef CAFE_COMMANDS_H_B2201016_8693_4379_9225_8FFEC3562AFE
#define CAFE_COMMANDS_H_B2201016_8693_4379_9225_8FFEC3562AFE

#include <vector>
#include <string>

int cafe_cmd_source(std::vector<std::string> tokens);

int cafe_shell_dispatch_command(char* cmd);

#endif

