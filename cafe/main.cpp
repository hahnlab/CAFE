extern "C" {
#include "cafe_shell.h"
#include <utils_string.h>
}

#include <time.h>

const char* __date__ = __DATE__;

const int VERBOSE = 0;

int main(int argc, char* argv[])
{
	//int i = 0;
	char prompt[STRING_BUF_SIZE];
	cafe_shell_init(VERBOSE);
	srand((unsigned int)time(NULL));
	int shell = 0;
	if ( argc == 2 )
	{
		cafe_cmd_source(argc,argv);
		cafe_cmd_exit(0,NULL);
	}
	else
	{
		while(shell !=  CAFE_SHELL_EXIT )
		{
			printf("%s", "# "); 
			if (!fgets(prompt,STRING_BUF_SIZE,stdin))
				return -1;
			shell = cafe_shell_dispatch_command(prompt);
		}
	}
  return 0;
}
