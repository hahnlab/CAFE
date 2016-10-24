#include <stdio.h>
#include <time.h>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <stdexcept>

#include "lambda.h"
#include "cafe_commands.h"
#include "reports.h"

extern "C" {
#include <utils_string.h>
#include "cafe_shell.h"
#include "cafe.h"

	extern pCafeParam cafe_param;
	extern void cafe_log(pCafeParam param, const char* msg, ...);
	void cafe_shell_clear_param(pCafeParam param, int btree_skip);

	extern pTree tmp_lambda_tree;
	extern pArrayList cafe_pCD;
	void __cafe_tree_string_gainloss(pString pstr, pPhylogenyNode ptnode);
	void __cafe_tree_string_sum_gainloss(pString pstr, pPhylogenyNode ptnode);
	pString __cafe_tree_gainloss_metapost(pCafeTree pcafe, int id, const char* title, double width, double height);
}

using namespace std;

typedef int(*cafe_command2)(pCafeParam cafe_param, vector<string>);

map<string, cafe_command2> get_dispatcher()
{
	map<string, cafe_command2> dispatcher;
	dispatcher["gainloss"] = cafe_cmd_gainloss;
	dispatcher["source"] = cafe_cmd_source;
	dispatcher["lambda"] = cafe_cmd_lambda;
	dispatcher["?"] = cafe_cmd_list;
	dispatcher["date"] = cafe_cmd_date;
	dispatcher["echo"] = cafe_cmd_echo;
	dispatcher["report"] = cafe_cmd_report;
	dispatcher["quit"] = cafe_cmd_exit;
	dispatcher["exit"] = cafe_cmd_exit;


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
		cafe_log(param, " %s", tokens[i].c_str());
	}
	cafe_log(param, "\n");
	return 0;
}

int cafe_cmd_date(pCafeParam param, vector<string> tokens)
{
	time_t now = time(NULL);
	struct tm *tm = localtime(&now);
	char buf[64];
	strftime(buf, sizeof(buf), "%a %b %e %T %Y", tm);
	cafe_log(param, "%s", buf);
	return 0;
}

int cafe_cmd_exit(pCafeParam param, vector<string> tokens)
{
	if (param->str_log)
	{
		string_free(param->str_log);
		fclose(param->flog);
		param->str_log = NULL;
	}
	if (tmp_lambda_tree) phylogeny_free(tmp_lambda_tree);
	tmp_lambda_tree = NULL;
	cafe_shell_clear_param(param, 0);
	if (cafe_pCD) arraylist_free(cafe_pCD, free);
	memory_free(param);
	return CAFE_SHELL_EXIT;
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

int cafe_cmd_gainloss(pCafeParam param, vector<string> tokens)
{
	if (param->pfamily == NULL)
		throw runtime_error("ERROR(gainloss): You did not load family: command 'load'\n");
	if (param->pcafe == NULL)
		throw runtime_error("ERROR(gainloss): You did not specify tree: command 'tree'\n");
	if (param->lambda == NULL)
		throw runtime_error("ERROR(gainloss): You did not set the parameters: command 'lambda' or 'lambdamu'\n");

	if (param->viterbi.viterbiNodeFamilysizes == NULL)
	{
		cafe_pCD = cafe_viterbi(param, cafe_pCD);
	}

	string name = tokens[1] + ".gs";
	FILE* fout = fopen(name.c_str(), "w");

	int i, j;
	int** nodefs = param->viterbi.viterbiNodeFamilysizes;
	pCafeTree pcafe = param->pcafe;
	int fsize = param->pfamily->flist->size;
	int nnodes = (pcafe->super.nlist->size - 1) / 2;

	pCafeTree psum = cafe_tree_copy(pcafe);

	for (i = 0; i < psum->super.nlist->size; i++)
	{
		pCafeNode pcnode = (pCafeNode)psum->super.nlist->array[i];
		pcnode->familysize = 0;
	}

	pString pstr;
	int sum = 0;
	int totalsum = 0;
	for (j = 0; j < psum->super.nlist->size; j++)
	{
		pCafeNode pcnode = (pCafeNode)psum->super.nlist->array[j];
		pcnode->viterbi[0] = 0;
		pcnode->viterbi[1] = 0;
	}

	for (i = 0; i < fsize; i++)
	{
		sum = 0;
		pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
		fprintf(fout, "%s\t", pitem->id);
		cafe_family_set_size(param->pfamily, i, pcafe);
		for (j = 0; j < nnodes; j++)
		{
			pCafeNode pnode = (pCafeNode)pcafe->super.nlist->array[2 * j + 1];
			pnode->familysize = nodefs[j][i];
		}
		for (j = 0; j < pcafe->super.nlist->size; j++)
		{
			pCafeNode pcnode = (pCafeNode)pcafe->super.nlist->array[j];
			pCafeNode pcsum = (pCafeNode)psum->super.nlist->array[j];
			if (tree_is_root((pTree)pcafe, (pTreeNode)pcnode)) continue;
			pCafeNode parent = (pCafeNode)((pTreeNode)pcnode)->parent;
			int diff = pcnode->familysize - parent->familysize;
			sum += diff;
			pcsum->familysize += diff;
			if (diff > 0)
			{
				pcsum->viterbi[0] += diff;
			}
			else if (diff <  0)
			{
				pcsum->viterbi[1] += diff;
			}
		}
		totalsum += sum;
		fprintf(fout, "%d\t", sum);
		pstr = phylogeny_string((pTree)pcafe, __cafe_tree_string_gainloss);
		fprintf(fout, "%s\n", pstr->buf);
		string_free(pstr);
	}
	fprintf(fout, "SUM\t%d\t", totalsum);
	pstr = phylogeny_string((pTree)psum, __cafe_tree_string_sum_gainloss);
	fprintf(fout, "%s\n", pstr->buf);
	string_free(pstr);

	fclose(fout);

	name = tokens[1] + ".mp";
	fout = fopen(name.c_str(), "w");

	ostringstream title;
	title << "SUM: " << totalsum;
	pstr = __cafe_tree_gainloss_metapost(psum, 0, title.str().c_str(), 10, 8);
	fprintf(fout, "%s", pstr->buf);
	string_free(pstr);
	cafe_tree_free(psum);
	fclose(fout);
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

	try
	{
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
	catch (std::exception& ex)
	{
		fprintf(stderr, "%s\n", ex.what());
		return -1;
	}
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

