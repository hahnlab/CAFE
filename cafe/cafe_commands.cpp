#include<stdio.h>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <iterator>
#include <algorithm>

#include "lambda.h"
#include "cafe_commands.h"

extern "C" {
#include <utils_string.h>
#include "cafe_shell.h"

	extern pCafeParam cafe_param;
	extern void cafe_log(pCafeParam param, char* msg, ...);
}

using namespace std;

typedef int(*cafe_command2)(pCafeParam cafe_param, vector<string>);

map<string, cafe_command2> get_dispatcher()
{
	map<string, cafe_command2> dispatcher;
	dispatcher["source"] = cafe_cmd_source;
	dispatcher["lambda"] = cafe_cmd_lambda;
	dispatcher["?"] = cafe_cmd_list;
	dispatcher["date"] = cafe_cmd_date;
	dispatcher["echo"] = cafe_cmd_echo;

	return dispatcher;
}

vector<string> tokenize(string s)
{
	vector<string> result;
	istringstream iss(s);

	while (iss.good()) {
		string tmp;
		iss >> tmp;
		if (tmp.size() > 0)
			result.push_back(tmp);
	}

	return result;
}

int cafe_cmd_echo(pCafeParam param, vector<string> tokens)
{
	for (size_t i = 1; i < tokens.size(); i++)
	{
		cafe_log(cafe_param, " %s", tokens[i].c_str());
	}
	cafe_log(cafe_param, "\n");
	return 0;
}

int cafe_cmd_date(pCafeParam param, vector<string> tokens)
{
	cafe_log(cafe_param, "%s", get_current_time());
	return 0;
}

int cafe_cmd_source(pCafeParam param, vector<string> tokens)
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

void list_commands(std::ostream& ost)
{
	vector<string> commands;

	int i;
	for (i = 0; cafe_cmd[i].command; i++)
	{
		commands.push_back(cafe_cmd[i].command);
	}
	map<string, cafe_command2> d = get_dispatcher();
	for (std::map<string, cafe_command2>::iterator iter = d.begin(); iter != d.end(); ++iter)
	{
		commands.push_back(iter->first);
	}
	sort(commands.begin(), commands.end());
	copy(commands.begin(), commands.end(), std::ostream_iterator<string>(ost, "\n"));
}

int cafe_cmd_list(pCafeParam, vector<string> tokens)
{
	list_commands(std::cout);
	return 0;
}



int cafe_shell_dispatch_command(char* cmd)
{
	using namespace std;

	map<string, cafe_command2> dispatcher = get_dispatcher();

	vector<string> tokens = tokenize(cmd);

	int i;
	int rtn = 0;

	if (tokens.size() > 0 && tokens[0][0] == '#')
		return 0;

	if (tokens.size() != 0)
	{
		rtn = CAFE_SHELL_NO_COMMAND;
		if (dispatcher.find(tokens[0]) != dispatcher.end())
			rtn = dispatcher[tokens[0]](cafe_param, tokens);
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

