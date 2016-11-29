#include <stdio.h>
#include <sys/stat.h>
#include <time.h>
#include <map>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <stdexcept>

#include "lambda.h"
#include "cafe_commands.h"
#include "reports.h"

/**
	\defgroup Commands Commands that are available in CAFE
*/

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
	void __cafe_cmd_viterbi_family_print(int idx);

	extern pBirthDeathCacheArray probability_cache;
}

using namespace std;

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
	dispatcher["genfamily"] = cafe_cmd_generate_random_family;
	dispatcher["log"] = cafe_cmd_log;
	dispatcher["version"] = cafe_cmd_version;
	dispatcher["info"] = cafe_cmd_print_param;
	dispatcher["load"] = cafe_cmd_load;
	dispatcher["family"] = cafe_cmd_family;

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

/**
\ingroup Commands
\brief Echoes test to the log file
*
*/
int cafe_cmd_echo(pCafeParam param, std::vector<std::string> tokens)
{
	for (size_t i = 1; i < tokens.size(); i++)
	{
		cafe_log(param, " %s", tokens[i].c_str());
	}
	cafe_log(param, "\n");
	return 0;
}

/**
\ingroup Commands
\brief Writes the current date and time to log file
*
*/
int cafe_cmd_date(pCafeParam param, std::vector<std::string> tokens)
{
	time_t now = time(NULL);
	struct tm *tm = localtime(&now);
	char buf[64];
	strftime(buf, sizeof(buf), "%a %b %e %T %Y", tm);
	cafe_log(param, "%s", buf);
	return 0;
}

/**
\ingroup Commands
\brief Close files, release memory and exit application
*
*/
int cafe_cmd_exit(pCafeParam param, std::vector<std::string> tokens)
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

/**
\ingroup Commands
\brief Sets file to which data is logged
*
* With no arguments, writes the current log file to stdout
* with the argument "stdout" the current log file is closed and log data goes to stdout
*
*/
int cafe_cmd_log(pCafeParam param, std::vector<std::string> tokens)
{
	if (tokens.size() == 1)
	{
		printf("Log: %s\n", param->flog == stdout ? "stdout" : param->str_log->buf);
	}
	else
	{
		string file_name;
		if (tokens[1] == "stdout")
			file_name = "stdout";
		else {
			ostringstream ost;
			copy(tokens.begin()+1, tokens.end(), std::ostream_iterator<string>(ost, " "));
			file_name = ost.str().substr(0, ost.str().size() - 1);
		}
		return set_log_file(param, file_name.c_str());
	}
	return 0;
}

void write_version(ostream &ost)
{
	ost << "Version: " << CAFE_VERSION << ", built at " << __DATE__ << "\n";
}

/**
\ingroup Commands
\brief Prints CAFE version and date of build
*
*/
int cafe_cmd_version(pCafeParam param, std::vector<std::string> tokens)
{
	write_version(cout);
	return 0;
}

/**
\ingroup Commands
\brief Logs various pieces of information about the application state
*
*/
int cafe_cmd_print_param(pCafeParam param, std::vector<std::string> tokens)
{
	log_param_values(param);
	return 0;
}


/**
\ingroup Commands
\brief Executes a series of commands from a CAFE command file
*
*/
int cafe_cmd_source(pCafeParam param, std::vector<std::string> tokens)
{
	if ( tokens.size() != 2 )
	{
		throw std::runtime_error("Usage: source <file>\n");
	}
	FILE* fp = fopen( tokens[1].c_str(), "r" );
	if ( fp == NULL ) 
	{
		ostringstream ost;
		ost << "Error(source): Cannot open " << tokens[1] << "\n";
		throw std::runtime_error(ost.str().c_str());
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

/**
\ingroup Commands
\brief List all commands available in CAFE
*
*/
int cafe_cmd_list(pCafeParam, std::vector<std::string> tokens)
{
	list_commands(std::cout);
	return 0;
}

void clear_node_viterbis(pTree ptree, pTreeNode ptnode, va_list ap1)
{
	pCafeNode pcnode = (pCafeNode)ptnode;
	pcnode->familysize = 0;
	pcnode->viterbi[0] = 0;
	pcnode->viterbi[1] = 0;

}

void clear_tree_viterbis(pCafeTree ptree)
{
	tree_traveral_infix((pTree)ptree, clear_node_viterbis);
}

void set_node_familysize(pCafeTree tree, int** node_family_sizes, int i)
{
	int nnodes = (tree->super.nlist->size - 1) / 2;
	for (int j = 0; j < nnodes; j++)
	{
		pCafeNode pnode = (pCafeNode)tree->super.nlist->array[2 * j + 1];
		pnode->familysize = node_family_sizes[j][i];
	}
}

int write_family_gainloss(ostream& ofst, std::string family_id, pCafeTree tree1, pCafeTree tree2)
{
	int sum = 0;
	ofst << family_id << "\t";
	for (int j = 0; j < tree1->super.nlist->size; j++)
	{
		pCafeNode pcnode = (pCafeNode)tree1->super.nlist->array[j];
		if (tree_is_root((pTree)tree1, (pTreeNode)pcnode)) continue;
		pCafeNode parent = (pCafeNode)((pTreeNode)pcnode)->parent;
		int diff = pcnode->familysize - parent->familysize;
		sum += diff;

		pCafeNode pcsum = (pCafeNode)tree2->super.nlist->array[j];
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

	ofst << sum << "\t";
	pString pstr = phylogeny_string((pTree)tree1, __cafe_tree_string_gainloss);
	ofst << pstr->buf << "\n";
	string_free(pstr);

	return sum;
}

const int REQUIRES_FAMILY = 0x01;
const int REQUIRES_TREE = 0x02;
const int REQUIRES_LAMBDA = 0x04;

void prereqs(pCafeParam param, int flags)
{
	if ((flags & REQUIRES_FAMILY) && param->pfamily == NULL)
		throw runtime_error("ERROR: You did not load family: command 'load'\n");
	if ((flags & REQUIRES_TREE) && param->pcafe == NULL)
		throw runtime_error("ERROR: You did not specify tree: command 'tree'\n");
	if ((flags & REQUIRES_LAMBDA) && param->lambda == NULL)
		throw runtime_error("ERROR: You did not set the parameters: command 'lambda' or 'lambdamu'\n");
}

/**
\ingroup Commands
\brief Write gains and losses
*
*/
int cafe_cmd_gainloss(pCafeParam param, std::vector<std::string> tokens)
{
	prereqs(param, REQUIRES_FAMILY | REQUIRES_TREE | REQUIRES_LAMBDA);

	if (param->viterbi.viterbiNodeFamilysizes == NULL)
	{
		cafe_pCD = cafe_viterbi(param, cafe_pCD);
	}

	string name = tokens[1] + ".gs";
	ofstream ofst(name.c_str());
	//FILE* fout = fopen(name.c_str(), "w");

	pCafeTree pcafe = param->pcafe;
	pCafeTree psum = cafe_tree_copy(pcafe);

	clear_tree_viterbis(psum);

	int totalsum = 0;
	int** nodefs = param->viterbi.viterbiNodeFamilysizes;
	int fsize = param->pfamily->flist->size;

	for (int i = 0; i < fsize; i++)
	{
		cafe_family_set_size(param->pfamily, i, pcafe);
		set_node_familysize(pcafe, nodefs, i);
		pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
		totalsum += write_family_gainloss(ofst, pitem->id, pcafe, psum);
	}
	ofst << "SUM\t" << totalsum << "\t";
	pString pstr = phylogeny_string((pTree)psum, __cafe_tree_string_sum_gainloss);
	ofst << pstr->buf << "\n";
	string_free(pstr);

	cafe_tree_free(psum);

	return 0;
}

vector<Argument> build_argument_list(vector<string> tokens)
{
	size_t i, j;
	vector<Argument> result;
	for (i = 1; i < tokens.size(); i++)
	{
		if (tokens[i][0] == '-' && !isdigit(tokens[i][1]))
		{
			Argument arg;
			arg.argc = 0;
			arg.opt = strdup(tokens[i].c_str());
			for (j = i + 1; j < tokens.size(); j++)
			{
				if (tokens[j][0] == '-' && !isdigit(tokens[j][1])) break;
				arg.argc++;
			}
			char ** argv = (char **)memory_new(tokens.size(), sizeof(char *));
			for (size_t k = 0, kk = i + 1; kk < tokens.size(); ++kk, ++k)
				argv[k] = strdup(tokens[kk].c_str());
			arg.argv = argv;
			result.push_back(arg);
			i = j - 1;
		}
	}
	return result;
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

int get_num_trials(vector<string> args)
{
	vector<string>::iterator it = find(args.begin(), args.end(), "-t");
	if (it == args.end())
		return 1;
	return atoi((++it)->c_str());
}

void verify_directory(string dir)
{
	if (dir.empty()) {    // check directory exists
		char * dirprefix = strdup(dir.c_str());
		char * pch = NULL;
		char * prevpch = NULL;
		char directory[STRING_BUF_SIZE];
		directory[0] = '\0';
		pch = strtok(dirprefix, "/\0");
		while (pch != NULL)
		{
			if (prevpch) {
				strcat(directory, prevpch);
				strcat(directory, "/");
			}
			prevpch = pch;
			pch = strtok(NULL, "/");
		}
		if (strlen(directory) != 0) {
			struct stat st;
			if (stat(directory, &st) != 0) {
				perror(directory);
				ostringstream ost;
				ost << "Please create directory " << dir << " before running genfamily.\n";
				throw std::runtime_error(ost.str().c_str());
			}
		}
	}
}

int* get_root_dist(pCafeTree pcafe, pCafeFamily pfamily, int k_value, int* family_sizes, int* rootfamily_sizes)
{
	int *root_dist = (int*)memory_new(pcafe->rfsize + 1, sizeof(int));
	int num_families = pfamily->flist->size;
	pCafeNode croot = (pCafeNode)pcafe->super.root;
	reset_birthdeath_cache(pcafe, k_value, family_sizes,rootfamily_sizes);
	printf("Viterbi\n");

	for (int i = 0; i < num_families; i++)
	{
		if (i % 1000 == 0)
		{
			printf("%d ...\n", i);
		}
		cafe_family_set_size(pfamily, i, pcafe);
		if (k_value > 0) {
			cafe_tree_clustered_viterbi(pcafe, k_value);
		}
		else {
			cafe_tree_viterbi(pcafe);
		}
		root_dist[croot->familysize]++;
	}

	return root_dist;
}

vector<int> get_clusters(int parameterized_k_value, int num_families, double* k_weights)
{
	int idx = 0;
	vector<int> result;
	if (parameterized_k_value > 0) {
		// assign clusters based on proportion (k_weights)
		int k_i;
		for (k_i = 0; k_i<parameterized_k_value - 1; k_i++) {
			for (int j = 0; j<k_weights[k_i] * num_families; j++) {
				result.push_back(k_i);
				idx++;
			}
		}
		for (; idx<num_families; idx++) {
			result.push_back(k_i);
		}
		// shuffle clusters
		random_shuffle(result.begin(), result.end());
	}

	return result;
}

void write_node_headers(ostream& s1, ostream& s2, pCafeTree pcafe)
{
	s1 << "DESC\tFID";
	s2 << "DESC\tFID";
	pArrayList nlist = pcafe->super.nlist;
	for (int i = 0; i < nlist->size; i += 2)
	{
		pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
		s1 << "\t" << pnode->name;
	}
	s1 << "\n";
	for (int i = 0; i < nlist->size; i++)
	{
		pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
		if (pnode->name) {
			s2 << "\t" << pnode->name;
		}
		else {
			s2 << "\t-" << i;
		}
	}
	s2 << "\n";
}

void write_leaves(ostream& ofst, pCafeTree pcafe, int *k, int i, int id, bool evens)
{
	pArrayList nlist = pcafe->super.nlist;

	// print test leaves
	if (k) {
		ofst << "k" << *k << "_";
	}
	ofst << "root" << i << "\t" << id;

	for (int n = 0; n < nlist->size; n += (evens ? 2 : 1))
	{
		pCafeNode pnode = (pCafeNode)nlist->array[n];
		ofst << "\t" << pnode->familysize;
	}
	ofst << "\n";
}

/**
\ingroup Commands
\brief Generates random families
*
*/
int cafe_cmd_generate_random_family(pCafeParam param, std::vector<std::string> tokens)
{
	if (tokens.size() == 1)
	{
		ostringstream ost;
		ost << "Usage: " << tokens[0] << " file\n";
		throw std::runtime_error(ost.str().c_str());
	}

	verify_directory(tokens[1]);

	int num_trials = get_num_trials(tokens);

	pCafeTree pcafe = param->pcafe;
	int j, n;

	prereqs(param, REQUIRES_TREE);

	int num_families = 0;
	if (param->root_dist == NULL) {
		prereqs(param, REQUIRES_FAMILY);

		num_families = param->pfamily->flist->size;
		param->param_set_func(param, param->parameters);
		param->root_dist = get_root_dist(pcafe, param->pfamily, param->parameterized_k_value, param->family_sizes, param->rootfamily_sizes);
	}
	else {
		num_families = 0;
		for (int i = 1; i <= pcafe->rfsize; i++)
		{
			num_families += param->root_dist[i];
		}
		prereqs(param, REQUIRES_LAMBDA);
		cafe_log(param, "Using user defined root size distribution for simulation... \n");
		param->param_set_func(param, param->parameters);
		reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, param->family_sizes, param->rootfamily_sizes);
	}

	int t;
	for (t = 0; t < num_trials; t++)
	{
		int id = 1;
		ostringstream ost;
		ost << tokens[1] << "_" << t + 1 << ".tab";
		ofstream ofst(ost.str().c_str());
		if (!ofst)
		{
			perror(tokens[1].c_str());
			ost << " failed to open";
			throw std::runtime_error(ost.str().c_str());
		}
		ostringstream tost;
		tost << tokens[1] << "_" << t + 1 << ".truth";
		ofstream truthst(tost.str().c_str());
		if (!truthst)
		{
			perror(tokens[1].c_str());
			tost << " failed to open";
			throw std::runtime_error(tost.str().c_str());
		}
		write_node_headers(ofst, truthst, pcafe);

		vector<int> k_arr = get_clusters(param->parameterized_k_value, num_families, param->k_weights);
		int idx = 0;
		for (int i = 1; i <= pcafe->rfsize; i++)
		{
			// iterates along root size distribution, and simulates as many families as the frequency of the root size.
			// need to make it pick random k based on k_weights. 
			if (param->root_dist[i])
			{
				pArrayList nlist = pcafe->super.nlist;
				for (j = 0; j < param->root_dist[i]; j++)
				{
					// cafe_tree_random_familysize uses birthdeath_matrix to calculate the probabilities.
					// if k > 0 point bd to k_bd[] before running cafe_tree_random_familysize to avoid EXC_BAD_ACCESS		
					if (param->parameterized_k_value > 0) {
						for (n = 0; n < nlist->size; n++)
						{
							pCafeNode pcnode = (pCafeNode)nlist->array[n];
							pcnode->birthdeath_matrix = birthdeath_cache_get_matrix(probability_cache, pcnode->super.branchlength, pcnode->birth_death_probabilities.lambda, pcnode->birth_death_probabilities.mu);
						}
					}
					// now randomly sample familysize
					cafe_tree_random_familysize(param->pcafe, i);

					int *k_value = param->parameterized_k_value > 0 ? &k_arr[idx] : NULL;
					write_leaves(ofst, pcafe, k_value, i, id, true);
					write_leaves(truthst, pcafe, k_value, i, id, false);

					id++;
					idx++;
				}
			}
		}
	}
	cafe_free_birthdeath_cache(pcafe);
	return 0;
}

struct load_args get_load_arguments(vector<Argument> pargs)
{
	struct load_args args;
	args.filter = false;
	args.num_threads = 0;
	args.num_random_samples = 0;
	args.pvalue = -1;

	for (size_t i = 0; i < pargs.size(); i++)
	{
		pArgument parg = &pargs[i];

		if (!strcmp(parg->opt, "-t"))
			sscanf(parg->argv[0], "%d", &args.num_threads);

		if (!strcmp(parg->opt, "-r"))
			sscanf(parg->argv[0], "%d", &args.num_random_samples);

		if (!strcmp(parg->opt, "-p"))
			sscanf(parg->argv[0], "%lf", &args.pvalue);

		if (!strcmp(parg->opt, "-l"))
		{
			if (!strcmp(parg->argv[0], "stdout"))
				args.log_file_name = "stdout";
			else {
				pString file_name = string_join(" ", parg->argc, parg->argv);
				args.log_file_name = file_name->buf;
				string_free(file_name);
			}
		}
		if ((!strcmp(parg->opt, "-filter")))
		{
			args.filter = true;
		}

		if ((!strcmp(parg->opt, "-i")))
		{
			pString file_name = string_join(" ", parg->argc, parg->argv);
			args.family_file_name = file_name->buf;
			string_free(file_name);
		}
	}

	return args;
}

void copy_args_to_param(pCafeParam param, struct load_args& args)
{
	if (args.num_threads > 0)
		param->num_threads = args.num_threads;
	if (args.num_random_samples > 0)
		param->num_random_samples = args.num_random_samples;
	if (args.pvalue > 0.0)
		param->pvalue = args.pvalue;
	if (!args.log_file_name.empty())
		set_log_file(param, args.log_file_name.c_str());
}

/**
\ingroup Commands
\brief Loads families from a family file with a defined format
*
* Takes six arguments: -t, -r, -p, -l, -i, and -filter
*/
int cafe_cmd_load(pCafeParam param, std::vector<std::string> tokens)
{
	if (tokens.size() < 2)
	{
		throw runtime_error("Usage(load): load <family file>\n");
	}
	cafe_shell_clear_param(param, 1);

	struct load_args args = get_load_arguments(build_argument_list(tokens));
	copy_args_to_param(param, args);

	if (args.filter && param->pcafe == NULL)
	{
		cerr << "Error(load): You did not specify tree. Skip filtering\n";
	}

	if (args.family_file_name.empty())
	{
		cafe_shell_clear_param(param, 1);
		throw runtime_error("ERROR(load): You must use -i option for input file\n");
	}

	param->str_fdata = string_new_with_string(args.family_file_name.c_str());
	param->pfamily = cafe_family_new(args.family_file_name.c_str(), 1);
	if (param->pfamily == NULL) 
		throw runtime_error("Failed to load file\n");

	// this is very important!!!!!!!!!!!
	param->rootfamily_sizes[0] = 1;
	//param->rootfamily_sizes[0] = 0;
	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	param->rootfamily_sizes[1] = MAX(30, rint(param->pfamily->max_size * 1.25));
	param->family_sizes[1] = param->pfamily->max_size + MAX(50, param->pfamily->max_size / 5);
	if (param->pcafe)
	{
		cafe_tree_set_parameters(param->pcafe, param->family_sizes, param->rootfamily_sizes, 0);
		cafe_family_set_species_index(param->pfamily, param->pcafe);
	}
	param->ML = (double*)memory_new(param->pfamily->flist->size, sizeof(double));
	param->MAP = (double*)memory_new(param->pfamily->flist->size, sizeof(double));

	if (param->pcafe && args.filter)
	{
		cafe_family_filter(param);
	}
	if (param->lambda && param->pcafe)
	{
		reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, param->family_sizes, param->rootfamily_sizes);
	}
	log_param_values(param);
	return 0;
}

struct family_args {
	int idx;
	string item_id;
	string add_id;
	vector<int> values;
	bool filter;
};

struct family_args get_family_arguments(vector<Argument> pargs)
{
	struct family_args args;
	args.idx = -1;
	args.filter = false;

	for (size_t i = 0; i < pargs.size(); i++)
	{
		pArgument parg = &pargs[i];

		if (!strcmp(parg->opt, "-idx"))
		{
			sscanf(parg->argv[0], "%d", &args.idx);
		}
		if (!strcmp(parg->opt, "-id"))
		{
			args.item_id = parg->argv[0];
		}
		if (!strcmp(parg->opt, "-filter"))
		{
			args.filter = true;
		}
		if (!strcmp(parg->opt, "-add"))
		{
			args.add_id = parg->argv[0];
			for (int j = 1; j <= parg->argc; j++)
			{
				int val;
				sscanf(parg->argv[j], "%d", &val);
				args.values.push_back(val);
			}
		}
	}

	return args;
}

int cafe_cmd_family(pCafeParam param, std::vector<std::string> tokens)
{
	prereqs(param, REQUIRES_FAMILY);

	int i, idx = 0;
	pCafeFamilyItem pitem = NULL;
	struct family_args args = get_family_arguments(build_argument_list(tokens));

	if (!args.item_id.empty())
	{
		args.idx = cafe_family_get_index(param->pfamily, args.item_id.c_str());
	}
	else if (!args.add_id.empty())
	{
		pCafeFamily pcf = param->pfamily;
		if (pcf == NULL)
		{
			pcf = (pCafeFamily)memory_new(1, sizeof(CafeFamily));
			pcf->flist = arraylist_new(11000);
			pArrayList nlist = param->pcafe->super.nlist;
			pcf->num_species = (nlist->size + 1) / 2;
			pcf->species = (char**)memory_new(pcf->num_species, sizeof(char*));
			pcf->index = (int*)memory_new(pcf->num_species, sizeof(int));
			for (i = 0; i < nlist->size; i += 2)
			{
				pcf->index[i] = i;
				pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
				pcf->species[i] = (char*)memory_new(strlen(pnode->name) + 1, sizeof(char));
				strcpy(pcf->species[i], pnode->name);
			}
			param->pfamily = pcf;
		}
		pCafeFamilyItem pitem = (pCafeFamilyItem)memory_new(1, sizeof(CafeFamilyItem));
		pitem->id = (char*)memory_new(args.add_id.length() + 1, sizeof(char));
		pitem->ref = -1;
		pitem->count = (int*)memory_new(pcf->num_species, sizeof(int));
		pitem->maxlh = -1;
		pitem->desc = NULL;
		strcpy(pitem->id, args.add_id.c_str());
		for (i = 1; i <= pcf->num_species; i++)
		{
			pitem->count[pcf->index[i]] = args.values[i];
		}
		arraylist_add(pcf->flist, pitem);
	}
	else if (args.filter)
	{
		prereqs(param, REQUIRES_TREE);

		cafe_family_filter(param);
		log_param_values(param);
		return 0;
	}
	if (idx < 0)
	{
		fprintf(stderr, "ERROR(family): your request not found\n");
		return -1;
	}
	if (idx >= param->pfamily->flist->size)
	{
		fprintf(stderr, "ERROR(family): The index range is from 0 to %d\n", param->pfamily->flist->size);
		return -1;
	}
	pitem = (pCafeFamilyItem)param->pfamily->flist->array[idx];
	if (pitem)
	{
		printf("ID: %s\n", pitem->id);
		printf("Desc: %s\n", pitem->desc);
		for (i = 0; i < param->pfamily->num_species; i++)
		{
			printf("%s: %d\n", param->pfamily->species[i], pitem->count[i]);
		}
		if (param->pcafe && probability_cache) __cafe_cmd_viterbi_family_print(idx);
	}

	return pitem ? 0 : -1;
}

