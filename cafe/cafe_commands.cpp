#include<stdio.h>
#include <map>
#include <string>

extern "C" {
#include <utils_string.h>
#include "cafe_shell.h"
}


int cafe_cmd_source(int argc, char* argv[])
{
	if ( argc != 2 )
	{
		fprintf( stderr, "Usage: %s <file>\n", argv[0]);
		return -1;
	}
	char fname[STRING_STEP_SIZE];
	string_pchar_join(fname, " ", argc-1, &argv[1]);
	FILE* fp = fopen( fname, "r" );
	if ( fp == NULL ) 
	{
		fprintf( stderr, "Error(source): Cannot open %s\n", argv[1] );
		return -1;
	}

	char buf[STRING_BUF_SIZE];
	int rtn  = 0;
	while( fgets( buf, STRING_BUF_SIZE, fp ) )
	{
		if ( (rtn = cafe_shell_dispatch_command(buf)) ) break;
	}
	fclose(fp);
	return rtn;
}


int cafe_shell_dispatch_command(char* cmd)
{
	using namespace std;

	map<string, cafe_command> dispatcher;
	dispatcher["source"] = cafe_cmd_source;
	pArrayList parg = string_pchar_space_split(cmd);

	int i;
	int rtn = 0;

	if (parg->size > 0)
	{
		if (((char*)parg->array[0])[0] == '#') return 0;
	}

	if (parg->size != 0)
	{
		rtn = CAFE_SHELL_NO_COMMAND;
		if (dispatcher.find((char *)parg->array[0]) != dispatcher.end())
			rtn = cafe_cmd[i].func(parg->size, (char**)parg->array);
		else for (i = 0; cafe_cmd[i].command; i++)
		{
			if (strcasecmp((char*)parg->array[0], cafe_cmd[i].command) == 0)
			{
				rtn = cafe_cmd[i].func(parg->size, (char**)parg->array);
				break;
			}
		}
		if (rtn == CAFE_SHELL_NO_COMMAND)
		{
			fprintf(stderr, "cafe: %s: command not found\n", (char*)parg->array[0]);
		}
	}
	arraylist_free(parg, NULL);
	return rtn;
}

extern "C" {
	int cafe_shell_dispatch_commandf(char* format, ...)
	{
		va_list ap;
		char buf[STRING_BUF_SIZE];
		va_start(ap, format);
		vsprintf(buf, format, ap);
		int r = cafe_shell_dispatch_command(buf);
		va_end(ap);
		return r;
	}
}

