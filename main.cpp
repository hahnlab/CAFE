#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <time.h>

extern "C" {
#include "cafe_shell.h"
#include <utils_string.h>

}
#include "cafe_commands.h"
#include "Globals.h"

const char* __date__ = __DATE__;

const int VERBOSE = 0;

int main(int argc, char* argv[])
{
	//int i = 0;
	char prompt[STRING_BUF_SIZE];
	Globals globals;
	srand((unsigned int)time(NULL));
	int shell = 0;
	if ( argc == 2 )
	{
		std::vector<std::string> tokens;
		tokens.push_back(argv[0]);
		tokens.push_back(argv[1]);
		cafe_cmd_source(globals, tokens);
		cafe_cmd_exit(globals, tokens);
	}
	else
	{
		while(shell !=  CAFE_SHELL_EXIT )
		{
			printf("%s", "# "); 
			if (!fgets(prompt,STRING_BUF_SIZE,stdin))
				return -1;
			shell = cafe_shell_dispatch_command(globals, prompt);
		}
	}
  return 0;
}
