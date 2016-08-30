#include "cafe_shell.h"

char* __date__ = __DATE__;

int main(int argc, char* argv[])
{
	//int i = 0;
	char prompt[STRING_BUF_SIZE];
	cafe_shell_init();
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
			printf("%s", "# "); fgets(prompt,STRING_BUF_SIZE,stdin);
			shell = cafe_shell_dispatch_command(prompt);
		}
	}
  return 0;
}
