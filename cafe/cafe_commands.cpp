#include<stdio.h>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>

extern "C" {
#include <utils_string.h>
#include "cafe_shell.h"
}

using namespace std;

vector<string> tokenize(string s)
{
	vector<string> result;
	istringstream iss(s);

	do
	{
		string sub;
		iss >> sub;
		result.push_back(sub);
	} while (iss);

	return result;
}

int cafe_cmd_source(vector<string> tokens)
{
	if ( tokens.size() != 2 )
	{
		fprintf( stderr, "Usage: source <file>\n");
		return -1;
	}
	FILE* fp = fopen( tokens[1].c_str(), "r" );
	if ( fp == NULL ) 
	{
		fprintf( stderr, "Error(source): Cannot open %s\n", tokens[1].c_str() );
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

typedef int(*cafe_command2)(vector<string>);

int cafe_shell_dispatch_command(char* cmd)
{
	using namespace std;

	map<string, cafe_command2> dispatcher;
	dispatcher["source"] = cafe_cmd_source;
	vector<string> tokens = tokenize(cmd);

	int i;
	int rtn = 0;

	if (tokens.size() > 0 && tokens[0][0] == '#')
		return 0;

	if (tokens.size() != 0)
	{
		rtn = CAFE_SHELL_NO_COMMAND;
		if (dispatcher.find(tokens[0]) != dispatcher.end())
			rtn = dispatcher[tokens[0]](tokens);
		else
		{
			pArrayList parg = string_pchar_space_split(cmd);
			for (i = 0; cafe_cmd[i].command; i++)
			{
				if (strcasecmp((char*)parg->array[0], cafe_cmd[i].command) == 0)
				{
					rtn = cafe_cmd[i].func(parg->size, (char**)parg->array);
					break;
				}
			}
			arraylist_free(parg, NULL);
		}
		if (rtn == CAFE_SHELL_NO_COMMAND)
		{
			fprintf(stderr, "cafe: %s: command not found\n", tokens[0].c_str());
		}
	}
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

