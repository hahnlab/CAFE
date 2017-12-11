#include "config.h"

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

#ifdef HAVE_READLINE_HISTORY_H
#include <readline/history.h>
#endif

#ifdef HAVE_READLINE_READLINE_H
#include <readline/readline.h>
#endif

const char* __date__ = __DATE__;

int main(int argc, char* argv[])
{
	//int i = 0;
	Globals globals;
	srand((unsigned int)time(NULL));
	int shell = 0;
	if ( argc == 2 )
	{
		std::vector<std::string> tokens;
		if (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0)
		{
			write_version(std::cout);
		}
		else
		{
			tokens.push_back(argv[0]);
			tokens.push_back(argv[1]);
			cafe_cmd_source(globals, tokens);
			cafe_cmd_exit(globals, tokens);
		}
	}
	else
	{
		while(shell !=  CAFE_SHELL_EXIT )
		{
#if defined(HAVE_READLINE_HISTORY_H) && defined(HAVE_READLINE_READLINE_H) 
			char *prompt = readline("# ");
			if (prompt == NULL)
				return -1;
			add_history(prompt); 
#else
			char prompt[STRING_BUF_SIZE];
			printf("%s", "# ");
			if (!fgets(prompt,STRING_BUF_SIZE,stdin))
				return -1;
#endif
			shell = cafe_shell_dispatch_command(globals, prompt);
		}
	}
  return 0;
}
