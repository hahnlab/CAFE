#include "../config.h"

#include <stdio.h>
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
#include <stack>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "lambda.h"
#include "lambdamu.h"
#include "cafe_commands.h"
#include "reports.h"
#include "pvalue.h"
#include "conditional_distribution.h"
#include "simerror.h"
#include "error_model.h"
#include "log_buffer.h"
#include "Globals.h"
#include "viterbi.h"
#include "gene_family.h"
#include "cross_validator.h"

#if HAVE_DIRENT_H == 1
#include<dirent.h>
#endif

#if HAVE_SYS_STAT_H == 1
#include <sys/stat.h>
#endif

/**
	\defgroup Commands Commands that are available in CAFE
*/

extern "C" {
#include <utils_string.h>
#include "cafe_shell.h"
#include "cafe.h"

	extern void cafe_log(pCafeParam param, const char* msg, ...);

	extern pTree tmp_lambda_tree;

	void __cafe_tree_string_gainloss(pString pstr, pPhylogenyNode ptnode);
	void __cafe_tree_string_sum_gainloss(pString pstr, pPhylogenyNode ptnode);

	extern pBirthDeathCacheArray probability_cache;
	int __cafe_cmd_extinct_count_zero(pTree pcafe);
	void __hg_print_sim_extinct(pHistogram** phist_sim_n, pHistogram* phist_sim,
		int r, pHistogram phist_tmp, double* cnt, int num_trials);
	pErrorStruct cafe_shell_create_error_matrix_from_estimate(pErrorMeasure errormeasure);
	void cafe_shell_set_branchlength(pCafeParam param, int max_family_size);
	void set_range_from_family(family_size_range* range, pCafeFamily family);
	double __cafe_best_lambda_search(double* plambda, void* args);
	double __cafe_cluster_lambda_search(double* parameters, void* args);
}

using namespace std;

/**
* \brief Holds the list of commands that are available in Cafe.
*
* Each element consists of a command and the function that is
* called to handle that command.
* Functions include #cafe_cmd_lambda, #cafe_cmd_family, #cafe_cmd_tree,
* etc.
*/
map<string, cafe_command2> get_dispatcher()
{
	map<string, cafe_command2> dispatcher;
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
	dispatcher["load"] = cafe_cmd_load;
	dispatcher["tree"] = cafe_cmd_tree;
	dispatcher["pvalue"] = cafe_cmd_pvalue;
	dispatcher["lhtest"] = cafe_cmd_lhtest;
	dispatcher["simerror"] = cafe_cmd_simerror;
	dispatcher["errormodel"] = cafe_cmd_errormodel;
	dispatcher["esterror"] = cafe_cmd_esterror;
	dispatcher["noerrormodel"] = cafe_cmd_noerrormodel;
	dispatcher["rootdist"] = cafe_cmd_rootdist;
	dispatcher["lambdamu"] = cafe_cmd_lambdamu;
    dispatcher["seed"] = cafe_cmd_seed;
#ifdef DEBUG
	dispatcher["score"] = cafe_cmd_score;
	dispatcher["simextinct"] = cafe_cmd_simextinct;
	dispatcher["branchlength"] = cafe_cmd_branchlength;
	dispatcher["accuracy"] = cafe_cmd_accuracy;
	dispatcher["gainloss"] = cafe_cmd_gainloss;
	dispatcher["info"] = cafe_cmd_print_param;
	dispatcher["cvfamily"] = cafe_cmd_cvfamily;
	dispatcher["cvspecies"] = cafe_cmd_cvspecies;
	dispatcher["extinct"] = cafe_cmd_extinct;
	dispatcher["family"] = cafe_cmd_family;
	dispatcher["retrieve"] = cafe_cmd_retrieve;
	dispatcher["save"] = cafe_cmd_save;
	dispatcher["viterbi"] = cafe_cmd_viterbi;
#endif


	return dispatcher;
}

void get_doubles_array(vector<double>& loc, pArgument parg)
{
	loc.resize(parg->argc);
	for (int j = 0; j < parg->argc; j++)
	{
		sscanf(parg->argv[j], "%lf", &loc[j]);
	}
}


io_error::io_error(string source, string file, bool write) : runtime_error("")
{
	ostringstream ost;
	ost << "ERROR (" << source << "): Cannot open " << file;
	if (write)
		ost << " in write mode.\n";
	text = ost.str();
}


/**
\ingroup Commands
\brief Echoes test to the log file
*
*/
int cafe_cmd_echo(Globals& globals, std::vector<std::string> tokens)
{
	log_buffer buf(&globals.param);
	ostream ost(&buf);
	for (size_t i = 1; i < tokens.size(); i++)
	{
		ost << " " << tokens[i].c_str();
	}
	ost << "\n";
	return 0;
}

/**
\ingroup Commands
\brief Writes the current date and time to log file
*
*/
int cafe_cmd_date(Globals& globals, std::vector<std::string> tokens)
{
	time_t now = time(NULL);
	struct tm *tm = localtime(&now);
	char buf[64];
	strftime(buf, sizeof(buf), "%a %b %e %T %Y", tm);
	cafe_log(&globals.param, "%s", buf);
	return 0;
}

/**
\ingroup Commands
\brief Close files, release memory and exit application
*
*/
int cafe_cmd_exit(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

	if (param->str_log)
	{
		string_free(param->str_log);
		fclose(param->flog);
		param->str_log = NULL;
	}
	if (tmp_lambda_tree) phylogeny_free(tmp_lambda_tree);
	tmp_lambda_tree = NULL;
	globals.Clear(0);

	return CAFE_SHELL_EXIT;
}

/**
\ingroup Commands
\ingroup Setters
\brief Sets file to which data is logged
*
* With no arguments, writes the current log file to stdout
* 
* with the argument "stdout" the current log file is closed and log data goes to stdout
*
*/
int cafe_cmd_log(Globals& globals, std::vector<std::string> tokens)
{
	if (tokens.size() == 1)
	{
		printf("Log: %s\n", globals.param.flog == stdout ? "stdout" : globals.param.str_log->buf);
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
		return set_log_file(&globals.param, file_name.c_str());
	}
	return 0;
}

void write_version(ostream &ost)
{
	ost << "Version: " << PACKAGE_VERSION << ", built at " << __DATE__ << "\n";
}

/**
\ingroup Commands
\brief Prints CAFE version and date of build
*
*/
int cafe_cmd_version(Globals& globals, std::vector<std::string> tokens)
{
	write_version(cout);
	return 0;
}

/**
\ingroup Commands
\brief Executes a series of commands from a CAFE command file
*
*/
int cafe_cmd_source(Globals& globals, std::vector<std::string> tokens)
{
	if ( tokens.size() != 2 )
	{
		throw std::runtime_error("Usage: source <file>\n");
	}

	FILE* fp = fopen( tokens[1].c_str(), "r" );
	if ( fp == NULL ) 
	{
		throw io_error("source", tokens[1], false);
	}

	char buf[STRING_BUF_SIZE];
	int rtn  = 0;
	while( fgets( buf, STRING_BUF_SIZE, fp ) )
	{
		if ( (rtn = cafe_shell_dispatch_command(globals, buf)) ) break;
	}
	fclose(fp);
	return rtn;
}

void list_commands(std::ostream& ost)
{
	vector<string> commands;

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
int cafe_cmd_list(Globals&, std::vector<std::string> tokens)
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

void prereqs(pCafeParam param, int flags)
{
	if ((flags & REQUIRES_FAMILY) && param->pfamily == NULL)
		throw runtime_error("ERROR: The gene families were not loaded. Please load gene families with the 'load' command.\n");
	if ((flags & REQUIRES_TREE) && param->pcafe == NULL)
		throw runtime_error("ERROR: The tree was not loaded. Please load a tree with the 'tree' command.\n");
	if ((flags & REQUIRES_LAMBDA) && param->lambda == NULL)
		throw runtime_error("ERROR: Lambda values were not set. Please set lambda values with the 'lambda' or 'lambdamu' commands.\n");
	if (flags & REQUIRES_ERRORMODEL)
	{
		pCafeFamily pcf = param->pfamily;
		if (!pcf->error_ptr && !pcf->errors)
			throw std::runtime_error("ERROR: No error model has been set. Please set an error model with the 'errormodel' command.\n");
	}
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
			arg.opt = tokens[i];
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


int cafe_shell_dispatch_command(Globals& globals, char* cmd)
{
	using namespace std;

	map<string, cafe_command2> dispatcher = get_dispatcher();

  vector<string> tokens = tokenize(cmd, REGULAR_WHITESPACE);

	int rtn = 0;

	if (tokens.size() > 0 && tokens[0][0] == '#')
		return 0;

	try
	{
		if (tokens.size() != 0)
		{
			rtn = CAFE_SHELL_NO_COMMAND;
			if (dispatcher.find(tokens[0]) != dispatcher.end())
				rtn = dispatcher[tokens[0]](globals, tokens);
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

int cafe_shell_dispatch_commandf(Globals& globals, const char* format, ...)
{
	va_list ap;
	char buf[STRING_BUF_SIZE];
	va_start(ap, format);
	vsprintf(buf, format, ap);
	int r = cafe_shell_dispatch_command(globals, buf);
	va_end(ap);
	return r;
}

int get_num_trials(vector<string> args)
{
	vector<string>::iterator it = find(args.begin(), args.end(), "-t");
	if (it == args.end())
		return 1;
	return atoi((++it)->c_str());
}

void cafe_shell_prompt(const char* prompt, const char* format, ...)
{
    va_list ap;
    va_start(ap, format);
    printf("%s ", prompt);
    if (vscanf(format, ap) == EOF)
        fprintf(stderr, "Read failure\n");
    va_end(ap);
}

int set_family_size_interactive(pCafeTree pcafe)
{
    int i;
    int max = 0;
    char buf[STRING_STEP_SIZE];

    pArrayList nlist = pcafe->super.nlist;
    for (i = 0; i < nlist->size; i += 2)
    {
        pCafeNode pnode = (pCafeNode)nlist->array[i];
        sprintf(buf, "%s: ", pnode->super.name);
        int size = -1;
        cafe_shell_prompt(buf, "%d", &size);
        if (size < 0)
        {
            fprintf(stderr, "ERROR: You put wrong data, you must enter an integer greater than or equal to 0\n");
            cafe_shell_prompt("Retry? [Y|N] ", "%s", buf);
            if (buf[0] != 'Y' && buf[0] != 'y') return -1;
            i -= 2;
        }
        else
        {
            pnode->familysize = size;
            if (size > max) max = size;
        }
    }
    return max;
}

void verify_directory(string dirname)
{
	if (dirname.empty())
		throw std::runtime_error("No directory name specified");

	size_t chr = dirname.find_last_of("/\\");
	if (chr == string::npos)
		throw std::runtime_error("Usage: genfamily directory/fileprefix -t integer");

	string directory = dirname.substr(0, chr);

#if HAVE_SYS_STAT_H == 1
	struct stat st;
	if (stat(directory.c_str(), &st) != 0) {
		perror(directory.c_str());
		ostringstream ost;
		ost << "Directory " << directory << " does not exist. It must be created before running genfamily.\n";
		throw std::runtime_error(ost.str().c_str());
	}
#endif

}

int* get_root_dist(pCafeTree pcafe, pCafeFamily pfamily, int k_value, family_size_range* range)
{
	int *root_dist = (int*)memory_new(pcafe->rfsize + 1, sizeof(int));
	int num_families = pfamily->flist->size;
	pCafeNode croot = (pCafeNode)pcafe->super.root;

	reset_birthdeath_cache(pcafe, k_value, range);
	printf("Viterbi\n");

	for (int i = 0; i < num_families; i++)
	{
		if (i % 1000 == 0)
		{
			printf("%d ...\n", i);
		}
    pCafeFamilyItem pitem = (pCafeFamilyItem)pfamily->flist->array[i];
		cafe_family_set_size(pfamily, pitem, pcafe);
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
int cafe_cmd_generate_random_family(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

	if (tokens.size() == 1)
	{
		throw std::runtime_error("Usage: genfamily directory/fileprefix -t integer");
	}

	verify_directory(tokens[1]);

	int num_trials = get_num_trials(tokens);

	pCafeTree pcafe = param->pcafe;
	int j, n;

	prereqs(param, REQUIRES_TREE);

	int num_families = 0;
	if (param->root_dist == NULL) {
		prereqs(param, REQUIRES_FAMILY | REQUIRES_LAMBDA);

		num_families = param->pfamily->flist->size;
		cafe_shell_set_lambdas(param, param->input.parameters);
		param->root_dist = get_root_dist(pcafe, param->pfamily, param->parameterized_k_value, &param->family_size);
	}
	else {
		num_families = 0;
		for (int i = 1; i <= pcafe->rfsize; i++)
		{
			num_families += param->root_dist[i];
		}
		prereqs(param, REQUIRES_LAMBDA);
		cafe_log(param, "Using user defined root size distribution for simulation... \n");
        cafe_shell_set_lambdas(param, param->input.parameters);
        reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, &param->family_size);
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
					cafe_tree_random_familysize(param->pcafe, i, probability_cache->maxFamilysize);

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
	args.max_size = -1;

	for (size_t i = 0; i < pargs.size(); i++)
	{
		pArgument parg = &pargs[i];

		if (!strcmp(parg->opt.c_str(), "-t"))
			sscanf(parg->argv[0], "%d", &args.num_threads);

		if (!strcmp(parg->opt.c_str(), "-r"))
			sscanf(parg->argv[0], "%d", &args.num_random_samples);

		if (!strcmp(parg->opt.c_str(), "-max_size"))
			sscanf(parg->argv[0], "%d", &args.max_size);

		if (!strcmp(parg->opt.c_str(), "-p"))
			sscanf(parg->argv[0], "%lf", &args.pvalue);

		if (!strcmp(parg->opt.c_str(), "-l"))
		{
			if (!strcmp(parg->argv[0], "stdout"))
				args.log_file_name = "stdout";
			else {
				pString file_name = string_join(" ", parg->argc, parg->argv);
				args.log_file_name = file_name->buf;
				string_free(file_name);
			}
		}
		if ((!strcmp(parg->opt.c_str(), "-filter")))
		{
			args.filter = true;
		}

		if ((!strcmp(parg->opt.c_str(), "-i")))
		{
			pString file_name = string_join(" ", parg->argc, parg->argv);
			args.family_file_name = file_name->buf;
			string_free(file_name);
		}
	}

	return args;
}

void copy_args_to_param(Globals& globals, struct load_args& args)
{
    if (args.num_threads > 0)
    {
        globals.param.num_threads = args.num_threads;
#if defined(_OPENMP)
        omp_set_num_threads(args.num_threads);
#endif
    }
	if (args.num_random_samples > 0)
		globals.num_random_samples = args.num_random_samples;
	if (args.pvalue > 0.0)
		globals.param.pvalue = args.pvalue;
	if (!args.log_file_name.empty())
		set_log_file(&globals.param, args.log_file_name.c_str());
}

bool endsWith(std::string const &fullString, std::string const &ending) {
  if (fullString.length() >= ending.length()) {
    return (0 == fullString.compare(fullString.length() - ending.length(), ending.length(), ending));
  }
  else {
    return false;
  }
}

/**
\ingroup Commands
\brief Loads families from a family file with a defined format
*
* Takes six arguments: -t, -r, -p, -l, -i, and -filter
*
* -filter ensures that there is at least one copy at the root (using parsimony) for each family. 
*/
int cafe_cmd_load(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

	if (tokens.size() < 2)
	{
		throw runtime_error("Usage(load): load <family file>\n");
	}
	globals.Clear(1);

	struct load_args args = get_load_arguments(build_argument_list(tokens));
	copy_args_to_param(globals, args);

	if (args.filter && param->pcafe == NULL)
	{
		cerr << "Error(load): You did not specify tree. Skip filtering\n";
	}

	if (args.family_file_name.empty())
	{
		globals.Clear(1);
		throw runtime_error("ERROR(load): You must use -i option for input file\n");
	}

  char separator = '\t';
  if (endsWith(args.family_file_name, "csv"))
    separator = ',';

	param->str_fdata = string_new_with_string(args.family_file_name.c_str());
  ifstream ifst(args.family_file_name);
	param->pfamily = load_gene_families(ifst, separator, args.max_size);
	if (param->pfamily == NULL) 
		throw runtime_error("Failed to load file\n");

	init_family_size(&param->family_size, param->pfamily->max_size);

	if (param->pcafe)
	{
		cafe_tree_set_parameters(param->pcafe, &param->family_size, 0);
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
		reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, &param->family_size);
	}
	log_buffer buf(&globals.param);
	ostream ost(&buf);
	log_param_values(ost, globals);
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

		if (!strcmp(parg->opt.c_str(), "-idx"))
		{
			sscanf(parg->argv[0], "%d", &args.idx);
		}
		if (!strcmp(parg->opt.c_str(), "-id"))
		{
			args.item_id = parg->argv[0];
		}
		if (!strcmp(parg->opt.c_str(), "-filter"))
		{
			args.filter = true;
		}
		if (!strcmp(parg->opt.c_str(), "-add"))
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

/**
\ingroup Commands
\brief Reports
*
*/
int cafe_cmd_report(Globals& globals, std::vector<std::string> tokens)
{
	prereqs(&globals.param, REQUIRES_FAMILY | REQUIRES_TREE | REQUIRES_LAMBDA);

	report_parameters params = get_report_parameters(tokens);

	cafe_do_report(globals, *globals.viterbi, &params);
	return 0;
}

template <typename T>
void write_vector(ostream& ost, vector<T> items, string delimiter)
{
	copy(items.begin(), items.end() - 1, ostream_iterator<T>(ost, delimiter.c_str()));
	ost << *(items.end() - 1);
}

void write_family(ostream& ost, pCafeFamily family)
{
	int i, n;
	
	vector<string> species(family->num_species);
	for (int i = 0; i < family->num_species; ++i)
		species[i] = family->species[i];
	
	ost << "Desc\tFamily ID\t";
	write_vector(ost, species, "\t");
	ost << "\n";

	for (i = 0; i < family->flist->size; i++)
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)family->flist->array[i];
		ost << pitem->desc << "\t" << pitem->id << "\t" << pitem->count[0];
		for (n = 1; n < family->num_species; n++)
		{
			ost << "\t" << pitem->count[n];
		}
		ost << "\n";
	}
}

/**
\ingroup Commands
\brief Tree
*
*/
int cafe_cmd_tree(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

	string newick;
	if (tokens.size() == 1)
	{
		printf("Newick: ");
		cin >> newick;
		if (cin.fail())
			fprintf(stderr, "Failed to read input\n");

	}
	else
	{
		ostringstream ost;
		copy(tokens.begin() + 1, tokens.end(), ostream_iterator<string>(ost));
		newick = ost.str();
	}
	if (param->pcafe)
	{
		if (probability_cache)
		{
			cafe_free_birthdeath_cache(param->pcafe);
		}
		cafe_tree_free(param->pcafe);
		memory_free(param->old_branchlength);
		param->old_branchlength = NULL;
	}
	param->pcafe = cafe_tree_new(newick.c_str(), &param->family_size, 0, 0);
	if (param->pcafe == NULL) {
		throw runtime_error("Failed to load tree from provided string");
	}
	param->num_branches = param->pcafe->super.nlist->size - 1;
	param->old_branchlength = (int*)memory_new(param->num_branches, sizeof(int));
	pTree ptree = (pTree)param->pcafe;
	// find max_branch_length and sum_branch_length.
	param->max_branch_length = 0;
	param->sum_branch_length = 0;
	for (int j = 0; j < ptree->nlist->size; j++)
	{
		pPhylogenyNode pnode = (pPhylogenyNode)ptree->nlist->array[j];
		if (pnode->branchlength > 0)
		{
			param->sum_branch_length += pnode->branchlength;
			if (param->max_branch_length < pnode->branchlength) {
				param->max_branch_length = pnode->branchlength;
			}
		}
		else if (!tree_is_root(ptree, (pTreeNode)pnode))
		{
			cafe_tree_free(param->pcafe);
			param->pcafe = NULL;
			throw runtime_error("Failed to load tree from provided string (branch length missing)");
		}
	}

	if (!param->quiet)
		cafe_tree_string_print(param->pcafe);
	if (param->pfamily)
	{
		cafe_family_set_species_index(param->pfamily, param->pcafe);
	}
	return 0;
}

struct viterbi_args get_viterbi_arguments(vector<Argument> pargs)
{
	struct viterbi_args args;
	args.idx = -1;
	args.all = false;

	for (size_t i = 0; i < pargs.size(); i++)
	{
		pArgument parg = &pargs[i];

		if (!strcmp(parg->opt.c_str(), "-all"))
		{
			args.all = true;
			if (parg->argc > 0)
			args.file = parg->argv[0];
		}
		if (!strcmp(parg->opt.c_str(), "-idx"))
		{
			sscanf(parg->argv[0], "%d", &args.idx);
			if (args.idx == -1)
			{
				throw std::runtime_error("ERROR(viterbi): idx parameter is not an integer\n");
			}
		}
		if (!strcmp(parg->opt.c_str(), "-id"))
		{
			args.item_id = parg->argv[0];
		}
	}

	return args;
}


int to_integer(std::string str)
{
	int size = -1;
	sscanf(str.c_str(), "%d", &size);
	if (size < 0)
		throw std::runtime_error("ERROR: Family sizes must be greater than or equal to 0\n");
	return size;
}

int cafe_shell_parse_familysize(pTree pcafe, std::vector<std::string> tokens)
{
	if (tokens.size() != (size_t)(pcafe->nlist->size / 2 + 1))
	{
		throw std::runtime_error("ERROR: Species count did not match value count\n");
	}
	vector<int> sizes(tokens.size());
	transform(tokens.begin(), tokens.end(), sizes.begin(), to_integer);
	for (size_t i = 0; i < sizes.size(); i++)
	{
		pCafeNode pnode = (pCafeNode)pcafe->nlist->array[i * 2];
		pnode->familysize = sizes[i];
	}
	return *max_element(sizes.begin(), sizes.end());
}

void run_viterbi_sim(pCafeTree pcafe, pCafeFamily pfamily, roots& roots)
{
	int familysize = pfamily->flist->size;
	roots.extinct.resize(familysize);
	roots.size.resize(familysize);
	roots.num.resize(pcafe->rfsize + 1);
	roots.avg_extinct.resize(pcafe->rfsize + 1);
	roots.total_extinct = 0;

	pCafeNode croot = (pCafeNode)pcafe->super.root;
	for (int i = 0; i < familysize; i++)
	{
    pCafeFamilyItem pitem = (pCafeFamilyItem)pfamily->flist->array[i];
    cafe_family_set_size(pfamily, pitem, pcafe);
		cafe_tree_viterbi(pcafe);
		roots.size[i] = croot->familysize;
		roots.extinct[i] = __cafe_cmd_extinct_count_zero((pTree)pcafe);
		roots.total_extinct += roots.extinct[i];
		roots.num[roots.size[i]]++;
	}

}

int init_histograms(int rfsize, roots& roots, int nsamples)
{
	roots.phist_data = (pHistogram*)memory_new(rfsize + 1, sizeof(pHistogram));
	roots.phist_sim = (pHistogram*)memory_new(rfsize + 1, sizeof(pHistogram));
	roots.phist_sim[0] = histogram_new(NULL, 0, 0);
	roots.phist_data[0] = histogram_new(NULL, 0, 0);
	int maxsize = roots.num[1];
	for (int i = 1; i <= rfsize; i++)
	{
		if (roots.num[i] > maxsize) maxsize = roots.num[i];
		if (roots.num[i])
		{
			roots.phist_data[i] = histogram_new(NULL, 0, 0);
			roots.phist_sim[i] = histogram_new(NULL, 0, 0);
		}
	}
	histogram_set_sparse_data(roots.phist_data[0], &roots.extinct.front(), nsamples);
	return maxsize;

}

ostream& operator<<(ostream& os, const Histogram& hist)
{
	int i;
	os << "MIN: " << hist.min << " ~ MAX: " << hist.max << "\n";
	os << "BIN: " << hist.nbins << "\n";
	os << "# Samples: " << hist.nsamples << "\n";
	if (hist.point)
	{
		for (i = 0; i < hist.nbins; i++)
		{
			os << hist.point[i] << "\t" << hist.count[i] << "\t" << ((double)hist.count[i]) / hist.nsamples << "\n";
		}
	}
	else
	{
		for (i = 0; i < hist.nbins; i++)
		{
			os << hist.min + hist.width*i << "\t" << hist.count[i] << "\t" << ((double)hist.count[i]) / hist.nsamples << "\n";
		}
	}

	return os;
}

struct pvalue_args get_pvalue_arguments(vector<Argument> pargs)
{
	struct pvalue_args args;
	args.index = -1;

	for (size_t i = 0; i < pargs.size(); i++)
	{
		pArgument parg = &pargs[i];

		if (!strcmp(parg->opt.c_str(), "-o"))
		{
			args.outfile = parg->argv[0];
		}
		if (!strcmp(parg->opt.c_str(), "-i"))
		{
			args.infile = parg->argv[0];
		}
		if (!strcmp(parg->opt.c_str(), "-idx"))
		{
			sscanf(parg->argv[0], "%d", &args.index);
			if (args.index == -1)
			{
				throw std::runtime_error("ERROR(pvalue): idx parameter is not an integer\n");
			}
		}
	}

	return args;
}

static pArrayList to_arraylist(matrix& v)
{
	pArrayList result = arraylist_new(10);
	for (size_t i = 0; i < v.size(); ++i)
	{
		double * temp = (double *)calloc(v[i].size(), sizeof(double));
		std::copy(temp, temp + v[i].size(), v[i].begin());
		arraylist_add(result, temp);
	}
	return result;
}


/**
\ingroup Commands
\brief Calculates pvalues

*/
int cafe_cmd_pvalue(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

	pvalue_args args = get_pvalue_arguments(build_argument_list(tokens));

	if (args.outfile.size() > 0)
	{
		prereqs(param, REQUIRES_TREE | REQUIRES_LAMBDA);
		ofstream ofst(args.outfile.c_str());
		if (!ofst)
		{
			throw io_error("pvalue", args.outfile, true);
		}
		matrix m = cafe_conditional_distribution(param->pcafe, &param->family_size, param->num_threads, globals.num_random_samples);
		pArrayList cond_dist = to_arraylist(m);
		write_pvalues(ofst, cond_dist, globals.num_random_samples);
		arraylist_free(cond_dist, free);
	}
	else if (args.infile.size() > 0)
	{
		cafe_log(param, "Loading p-values ... \n");
		ifstream ifst(args.infile.c_str());
		if (!ifst)
		{
			throw io_error("pvalue", args.outfile, true);
		}
		read_pvalues(ifst, globals.num_random_samples);
		cafe_log(param, "Done Loading p-values ... \n");
	}
	else if (args.index >= 0)
	{
		pvalues_for_family(param->pcafe, param->pfamily, &param->family_size, param->num_threads, globals.num_random_samples, args.index);
	}
	else if (tokens.size() > 1)
	{
		prereqs(param, REQUIRES_TREE);
		print_pvalues(cout, param->pcafe, cafe_shell_parse_familysize((pTree)param->pcafe, tokens), globals.num_random_samples, probability_cache);
	}
	else
	{
		prereqs(param, REQUIRES_TREE);
		print_pvalues(cout, param->pcafe, set_family_size_interactive(param->pcafe), globals.num_random_samples, probability_cache);
	}
	return 0;
}

lhtest_args get_lhtest_arguments(vector<Argument> pargs)
{
	lhtest_args args;
	args.lambda = 0.0;

	for (size_t i = 0; i < pargs.size(); i++)
	{
		pArgument parg = &pargs[i];

		if (!strcmp(parg->opt.c_str(), "-d"))
		{
			args.directory = parg->argv[0];
		}
		if (!strcmp(parg->opt.c_str(), "-l"))
		{
			sscanf(parg->argv[0], "%lf", &args.lambda);
		}
		if (!strcmp(parg->opt.c_str(), "-t"))
		{
			args.tree = parg->argv[0];
		}
		if (!strcmp(parg->opt.c_str(), "-o"))
		{
			args.outfile = parg->argv[0];
		}
	}

	return args;
}


std::vector<std::string> enumerate_directory(std::string dir, std::string ext)
{
  std::vector<std::string> pal;
#if HAVE_DIRENT_H == 1
  DIR* d = opendir(dir.c_str());
	struct dirent* dp;
	if (d == NULL) 
		throw std::runtime_error("Failed to read directory");

	while ((dp = readdir(d)) != NULL)
	{
		if (!ext.empty())
		{
			char* point = rindex(dp->d_name, '.');
			if (point == NULL) continue;
			point++;
			if (strcmp(point, ext.c_str())) continue;
		}
		pal.push_back(dp->d_name);
	}
	closedir(d);
#else
  throw std::runtime_error("Unable to read directory on this system");
#endif

	return pal;
}

int cafe_cmd_lhtest(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

	lhtest_args args = get_lhtest_arguments(build_argument_list(tokens));

	FILE* fout = stdout;
	if (!args.outfile.empty())
	{
		if ((fout = fopen(args.outfile.c_str(), "w")) == NULL)
		{
			perror(args.outfile.c_str());
			throw std::runtime_error("Failed to open output file");
		}
	}
	std::vector<std::string> files = enumerate_directory(args.directory, "tab");
	pString pstr_cafe = phylogeny_string((pTree)param->pcafe, NULL);
	int j;
    int num_families = param->pfamily->flist->size;
    std::vector<double> ml(num_families), map(num_families), pr(FAMILYSIZEMAX);
    copy(param->ML, param->ML + num_families, ml.begin());
    copy(param->MAP, param->MAP + num_families, map.begin());
    copy(param->prior_rfsize, param->prior_rfsize + FAMILYSIZEMAX, pr.begin());

	for (size_t i = 0; i < files.size(); i++)
	{
		std::string fname = files[i];
		if (fname[0] == '.') continue;
		cafe_shell_dispatch_commandf(globals, "load -i %s/%s -p 0.01 -t 10 -l %s",
			args.directory.c_str(), fname.c_str(), param->str_log ? param->str_log->buf : "stdout");
		cafe_shell_dispatch_commandf(globals, "tree %s", pstr_cafe->buf);
		cafe_shell_dispatch_commandf(globals, "lambda -s -l %lf", args.lambda);
		try
		{
			double p = get_posterior(param->pfamily, param->pcafe, &param->family_size, ml, map, pr, param->quiet);
			fprintf(fout, "\t%lf\t%lf", p, param->lambda[0]);
		}
		catch (std::runtime_error& e)
		{
			cerr << e.what() << endl;
			fprintf(fout, "\t%lf\t%lf", log(0), param->lambda[0]);
		}
		cafe_family_reset_maxlh(param->pfamily);
		cafe_shell_dispatch_commandf(globals, "lambda -s -v %lf -t %s", args.lambda, args.tree.c_str());
		try
		{
			double p = get_posterior(param->pfamily, param->pcafe, &param->family_size, ml, map, pr, param->quiet);
			fprintf(fout, "\t%lf", p);
		}
		catch (std::runtime_error& e)
		{
			cerr << e.what() << endl;
			fprintf(fout, "\t%lf", log(0));
		}

		for (j = 0; j < param->num_lambdas; j++)
		{
			fprintf(fout, "\t%lf", param->lambda[j]);
		}
		fprintf(fout, "\n");
		fflush(fout);
		globals.Clear(0);
	}
	fclose(fout);
	string_free(pstr_cafe);
	return 0;
}

int cafe_cmd_simerror(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

	prereqs(param, REQUIRES_FAMILY | REQUIRES_TREE | REQUIRES_LAMBDA | REQUIRES_ERRORMODEL);

	std::vector<Argument> args = build_argument_list(tokens);

	std::string prefix;
	int repeat = 1;
	for (size_t i = 0; i < args.size(); i++)
	{
		// Search for whole family 
		if (!strcmp(args[i].opt.c_str(), "-pre"))
		{
			for (int j = 0; j < args[i].argc; ++j)
				prefix += args[i].argv[j];
		}
		else if (!strcmp(args[i].opt.c_str(), "-rep"))
		{
			repeat = atoi(args[i].argv[0]);
		}
	}
	simerror(param->pfamily, prefix, repeat);

	return 0;
}

struct errormodel_args
{
	std::string model_file;
	std::vector<std::string> species_list;
	bool all;
};

errormodel_args get_errormodel_arguments(vector<Argument> pargs)
{
	errormodel_args args;
	args.all = false;

	for (size_t i = 0; i < pargs.size(); i++)
	{
		pArgument parg = &pargs[i];

		if (!strcmp(parg->opt.c_str(), "-model"))
		{
			args.model_file = parg->argv[0];
		}
		if (!strcmp(parg->opt.c_str(), "-sp"))
		{
			for (int j = 0; j < parg->argc; j++)
			{
				args.species_list.push_back(parg->argv[j]);
			}
		}
		if (!strcmp(parg->opt.c_str(), "-all"))
		{
			args.all = true;
		}
	}

	return args;
}

int cafe_cmd_errormodel(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

	prereqs(param, REQUIRES_FAMILY | REQUIRES_TREE);

	errormodel_args args = get_errormodel_arguments(build_argument_list(tokens));

	if (!args.model_file.empty())
	{
		if (!args.species_list.empty())
		{
			for (size_t j = 0; j < args.species_list.size(); ++j)
			{
				set_error_matrix_from_file(param->pfamily, param->pcafe, param->family_size, args.model_file, args.species_list[j]);
			}
		}
		else if (args.all)
		{
			set_error_matrix_from_file(param->pfamily, param->pcafe, param->family_size, args.model_file, std::string());
		}
		fprintf(stderr, "errormodel: %s set.\n", args.model_file.c_str());
		fprintf(stderr, "errormodel: Remember that the rows in the errormodel file have to add up to 1 (rows in the errormodel file correspond to columns in the errormatrix).\n");
		fprintf(stderr, "errormodel: The program does not check, only renormalizes.\n");
	}

	if (!param->pfamily->error_ptr || !param->pfamily->errors) {
		throw std::runtime_error("ERROR(errormodel): we need an error model specified (-model) or two data files.\n");
	}

	return 0;
}

esterror_args get_esterror_arguments(vector<Argument> pargs)
{
	esterror_args args;
	args.symmetric = false;
	args.peakzero = false;
	args.max_diff = 2;

	for (size_t i = 0; i < pargs.size(); i++)
	{
		pArgument parg = &pargs[i];

		if (!strcmp(parg->opt.c_str(), "-o"))
		{
			args.outfile = parg->argv[0];
		}
		if (!strcmp(parg->opt.c_str(), "-datatrue"))
		{
			args.truth_file = parg->argv[0];
		}
		if (!strcmp(parg->opt.c_str(), "-dataerror"))
		{
			for (int j = 0; j < parg->argc; j++)
			{
				args.data_error_files.push_back(parg->argv[j]);
			}
		}
		if (!strcmp(parg->opt.c_str(), "-symm"))
		{
			args.symmetric = true;
		}
		if (!strcmp(parg->opt.c_str(), "-diff"))
		{
			args.max_diff = atoi(parg->argv[0]);
		}
		if (!strcmp(parg->opt.c_str(), "-peakzero"))
		{
			args.peakzero = true;
		}
	}

	return args;
}

void validate(esterror_args args)
{
	if (args.outfile.empty())
		throw std::runtime_error("ERROR(esterror): need to specify the output file \"-o outfile\".\n");

	size_t num_error_files = args.data_error_files.size();
	if (num_error_files < 1 || num_error_files > 2)
		throw std::runtime_error("ERROR(esterror): [-dataerror file1 file2] or [-dataerror file1 -datatrue file2] -o outfile.\n");

	if (num_error_files == 1 && args.truth_file.empty())
	{
		throw std::runtime_error("ERROR(esterror): we need another data file with error or a true data file to compare.\n");
	}

}

int cafe_cmd_esterror(Globals& globals, std::vector<std::string> tokens)
{
	esterror_args args = get_esterror_arguments(build_argument_list(tokens));

	validate(args);

	ofstream ofst(args.outfile.c_str());
	if (!ofst)
	{
		throw io_error("esterror", args.outfile, true);
	}

	pErrorMeasure errormeasure = NULL;

	log_buffer buf(&globals.param);
	ostream ost(&buf);
	// estimate error matrix                
	if (args.data_error_files.size() == 2) {
        std::ifstream ist1(args.data_error_files[0]), ist2(args.data_error_files[1]);
        errormeasure = estimate_error_double_measure(ost, &ist1, &ist2, args.symmetric, args.max_diff, args.peakzero, globals.param.family_size.max);
	}
	else if (args.data_error_files.size() == 1) {
		errormeasure = estimate_error_true_measure(ost, args.data_error_files[0].c_str(), args.truth_file.c_str(), args.symmetric, args.max_diff, args.peakzero, globals.param.family_size.max);
	}

	// create errormodel based on errormeasure
	pErrorStruct errormodel = cafe_shell_create_error_matrix_from_estimate(errormeasure);

	// write errormodel
	ofst << *errormodel;

	return 0;

}

int cafe_cmd_noerrormodel(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

	prereqs(param, REQUIRES_TREE | REQUIRES_FAMILY);
	std::string species;
	bool all = false;

	vector<Argument> args = build_argument_list(tokens);
	for (size_t i = 0; i < args.size(); i++)
	{
		pArgument parg = &args[i];

		if (!strcmp(parg->opt.c_str(), "-sp"))
		{
			ostringstream ost;
			copy(parg->argv, parg->argv+ parg->argc, std::ostream_iterator<string>(ost, " "));
			species = ost.str();
		}
		if (!strcmp(parg->opt.c_str(), "-all"))
		{
			all = true;
		}
	}
	
	if ((tokens.size() >= 3) && !species.empty())
	{
		remove_error_model(param->pfamily, param->pcafe, species);
	}
	else if (tokens.size() == 2 && all)
	{
		free_error_model(param->pfamily, param->pcafe);
	}
	return 0;
}

void tree_set_branch_lengths(pCafeTree pcafe, std::vector<int> lengths)
{
	pArrayList nlist = pcafe->super.nlist;

	if (lengths.size() != (size_t)nlist->size)
	{
		std::ostringstream ost;
		ost << "ERROR: There are " << nlist->size << " branches including the empty branch of root\n";
		throw std::runtime_error(ost.str());
	}
	for (size_t i = 0; i < lengths.size(); i++)
	{
		pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
		if (lengths[i] > 0)
		{
			pnode->branchlength = lengths[i];
		}
		else
		{
			fprintf(stderr, "ERROR: the branch length of node %d is not changed\n", (int)i);
		}
	}
	if (probability_cache) cafe_tree_set_birthdeath(pcafe, probability_cache->maxFamilysize);

}

int s_to_i(std::string s)
{
	int size = -1;
	sscanf(s.c_str(), "%d", &size);
	return size;
}

std::string get_input_file(std::vector<std::string> tokens)
{
	vector<Argument> args = build_argument_list(tokens);
	for (size_t i = 0; i < args.size(); i++)
	{
		pArgument parg = &args[i];

		if (!strcmp(parg->opt.c_str(), "-i"))
		{
      return parg->argv[0];
		}
	}

	throw std::runtime_error("No input file specified");
}

/**
\ingroup Commands
\brief Specify root family size distribution for simulation

Arguments: -i input file.
*/
int cafe_cmd_rootdist(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

	prereqs(param, REQUIRES_TREE);

	std::string file = get_input_file(tokens);
	if (tokens.size() < 2)
	{
		prereqs(param, REQUIRES_FAMILY|REQUIRES_LAMBDA);

		cafe_log(param, "-----------------------------------------------------------\n");
		cafe_log(param, "Family information: %s\n", param->str_fdata->buf);
		cafe_log(param, "Log: %s\n", param->flog == stdout ? "stdout" : param->str_log->buf);
		if (param->pcafe)
		{
			pString pstr = phylogeny_string((pTree)param->pcafe, NULL);
			cafe_log(param, "Tree: %s\n", pstr->buf);
			string_free(pstr);
		}
		if (param->lambda)
		{
			pString pstr = cafe_tree_string_with_lambda(param->pcafe);
			cafe_log(param, "Lambda: %s\n", pstr->buf);
			string_free(pstr);
		}
		cafe_log(param, "The number of families is %d\n", param->pfamily->flist->size);
		int i;
		pCafeTree pcafe = param->pcafe;

		reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, &param->family_size);
		for (i = 0; i< param->pfamily->flist->size; i++)
		{
      pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
      cafe_family_set_size(param->pfamily, pitem, pcafe);
			cafe_tree_viterbi(pcafe);
			cafe_log(param, "%d\n", ((pCafeNode)pcafe->super.root)->familysize);
		}
		cafe_free_birthdeath_cache(pcafe);
		cafe_log(param, "\n");
	}
	else if (!file.empty())
	{
		FILE* fp = fopen(file.c_str(), "r");
    char buf[STRING_BUF_SIZE];
		if (fp == NULL)
		{
			throw io_error("rootdist", file, false);
		}
		if (fgets(buf, STRING_BUF_SIZE, fp) == NULL)
		{
			fclose(fp);
			fprintf(stderr, "Empty file: %s\n", file.c_str());
			return -1;
		}
		int i = 0;
		int max_rootsize = 0;
		string_pchar_chomp(buf);
		pArrayList data = string_pchar_split(buf, ' ');
		pArrayList max = string_pchar_split((char *)data->array[data->size - 1], ':');
		max_rootsize = atoi((char*)max->array[1]);
		arraylist_free(data, NULL);
		if (param->root_dist) { memory_free(param->root_dist); param->root_dist = NULL; }
		param->root_dist = (int*)memory_new(max_rootsize + 1, sizeof(int));

		param->family_size.root_min = 1;
		param->family_size.root_max = max_rootsize;
		param->family_size.min = 0;
		param->family_size.max = max_rootsize * 2;
		copy_range_to_tree(param->pcafe, &param->family_size);

		for (i = 0; fgets(buf, STRING_BUF_SIZE, fp); i++)
		{
			string_pchar_chomp(buf);
			data = string_pchar_split(buf, ' ');
			param->root_dist[atoi((char *)data->array[0])] = atoi((char *)data->array[1]);
		}
		arraylist_free(data, NULL);
	}

	return 0;
}

void log_param_values(std::ostream& ost, Globals& globals)
{
	pCafeParam param = &globals.param;

	ost << "-----------------------------------------------------------\n";
	ost << "Family information: " << param->str_fdata->buf << "\n";
	ost << "Log: " << (param->flog == stdout ? "stdout" : param->str_log->buf) << "\n";
	if (param->pcafe)
	{
		pString pstr = phylogeny_string((pTree)param->pcafe, NULL);
		ost << "Tree: " << pstr->buf << "\n";
		string_free(pstr);
	}
	ost << "The number of families is " << param->pfamily->flist->size << "\n";
	ost << "Root Family size : " << param->family_size.root_min << " ~ " << param->family_size.root_max << "\n";
	ost << "Family size : " << param->family_size.min << " ~ " << param->family_size.max << "\n";
	ost << "P-value: " << param->pvalue << "\n";
	ost << "Num of Threads: " << param->num_threads << "\n";
	ost << "Num of Random: " << globals.num_random_samples << "\n";
	if (param->lambda)
	{
		pString pstr = cafe_tree_string_with_lambda(param->pcafe);
		ost << "Lambda: " << pstr->buf << "\n";
		string_free(pstr);
	}
}

/**
\ingroup Commands
\brief Set the random seed for reproducible results
*
*/
int cafe_cmd_seed(Globals& globals, std::vector<std::string> tokens)
{
    if (tokens.size() < 2)
    {
        throw std::runtime_error("No value provided for seed");
    }
    try
    {
        srand(std::stoi(tokens[1]));
    }
    catch (...)
    {
        throw std::runtime_error("Failed to set seed from value " + tokens[1]);
    }
    return 0;
}

#ifdef DEBUG

void viterbi_print(pCafeTree pcafe, int max)
{
	check_cache_and_compute_likelihoods(pcafe, max, probability_cache);
	double* lh = get_likelihoods(pcafe);
	double mlh = __max(lh, pcafe->rfsize);
	cafe_tree_viterbi(pcafe);
	pString pstr = cafe_tree_string(pcafe);
	printf("%g\t%s\n", mlh, pstr->buf);
	string_free(pstr);
}

double cafe_shell_score(Globals& globals)
{
	int i = 0;
	double score = 0;
	if (globals.param.parameterized_k_value > 0) {
		if (globals.param.num_mus > 0) {
			score = -cafe_cluster_lambda_mu_search(globals.param.input.parameters, (void*)&globals.param);
			// print
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			for (i = 0; i<globals.param.num_lambdas; i++) {
				string_pchar_join_double(buf, ",", globals.param.parameterized_k_value, &globals.param.input.parameters[i*globals.param.parameterized_k_value]);
				cafe_log(&globals.param, "Lambda branch %d: %s\n", i, buf);
				buf[0] = '\0';
			}
			for (i = 0; i<globals.param.num_mus; i++) {
				string_pchar_join_double(buf, ",", globals.param.parameterized_k_value, &globals.param.input.parameters[globals.param.num_lambdas*globals.param.parameterized_k_value + i*globals.param.parameterized_k_value]);
				cafe_log(&globals.param, "Mu branch %d: %s \n", i, buf);
				buf[0] = '\0';
			}
			if (globals.param.parameterized_k_value > 0) {
				string_pchar_join_double(buf, ",", globals.param.parameterized_k_value, globals.param.k_weights);
				cafe_log(&globals.param, "p : %s\n", buf);
			}
			cafe_log(&globals.param, "Score: %f\n", score);

		}
		else {
			score = -__cafe_cluster_lambda_search(globals.param.input.parameters, (void*)&globals.param);
			// print
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			string_pchar_join_double(buf, ",", globals.param.num_lambdas*globals.param.parameterized_k_value, globals.param.input.parameters);
			cafe_log(&globals.param, "Lambda : %s\n", buf);
			buf[0] = '\0';
			if (globals.param.parameterized_k_value > 0) {
				string_pchar_join_double(buf, ",", globals.param.parameterized_k_value, globals.param.k_weights);
				cafe_log(&globals.param, "p : %s\n", buf);
			}
			cafe_log(&globals.param, "Score: %f\n", score);
		}
	}
	else {
		if (globals.param.num_mus > 0) {
			score = -cafe_best_lambda_mu_search(globals.param.input.parameters, (void*)&globals.param);
			// print
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			string_pchar_join_double(buf, ",", globals.param.num_lambdas, globals.param.input.parameters);
			cafe_log(&globals.param, "Lambda : %s ", buf, score);
			buf[0] = '\0';
			string_pchar_join_double(buf, ",", globals.param.num_mus, globals.param.input.parameters + globals.param.num_lambdas);
			cafe_log(&globals.param, "Mu : %s & Score: %f\n", buf, score);
		}
		else {
			score = -__cafe_best_lambda_search(globals.param.input.parameters, (void*)&globals.param);
			// print
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			string_pchar_join_double(buf, ",", globals.param.num_lambdas, globals.param.input.parameters);
			cafe_log(&globals.param, "Lambda : %s & Score: %f\n", buf, score);
		}
	}
	return score;
}

/**
\ingroup Commands
\brief Score
*
*/
int cafe_cmd_score(Globals& globals, std::vector<std::string> tokens)
{
  log_buffer buf(&globals.param);
  ostream ost(&buf);
  
  pCafeParam param = &globals.param;

	double score = cafe_shell_score(globals);
	ost << score << endl;
	if (param->parameterized_k_value > 0) {
		log_cluster_membership(param->pfamily, param->parameterized_k_value, param->p_z_membership, ost);
	}
	cafe_tree_set_parameters(param->pcafe, &param->family_size, 0);
	return 0;
}

/**
\ingroup Commands
\brief Runs a Monte Carlo simulation against the data and reports the number of extinctions that occurred.

Arguments: -r range Can be specified as a max or a colon-separated range
-t Number of trials to run
*/
int cafe_cmd_simextinct(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

	prereqs(param, REQUIRES_TREE | REQUIRES_LAMBDA);

	int num_trials = 10000;
	int range[2] = { 1, param->family_size.root_max };
	vector<Argument> args = build_argument_list(tokens);
	for (size_t i = 0; i < args.size(); i++)
	{
		pArgument parg = &args[i];

		if (!strcmp(parg->opt.c_str(), "-t"))
		{
			sscanf(parg->argv[0], "%d", &num_trials);
		}

		if (!strcmp(parg->opt.c_str(), "-r"))
		{
			if (std::string(parg->argv[0]).find(':') != std::string::npos)
			{
				sscanf(parg->argv[0], "%d:%d", &range[0], &range[1]);
			}
			else
			{
				sscanf(parg->argv[0], "%d", &range[1]);
				range[0] = range[1];
			}
		}
	}

	cafe_log(param, "Extinction count from Monte Carlo:\n");
	cafe_log(param, "root range: %d ~ %d\n", range[0], range[1]);
	cafe_log(param, "# trials: %d\n", num_trials);

	if (range[0] > range[1] || range[1] > param->family_size.root_max)
	{
		ostringstream ost;
		ost << "ERROR(simextinct): -r : 1 ~ " << param->family_size.root_max << "\n";
		throw std::runtime_error(ost.str());
	}

	int i, r;
	unsigned int accu_sum = 0;
	pHistogram phist_sim = histogram_new(NULL, 0, 0);
	pHistogram phist_accu = histogram_new(NULL, 0, 0);
	vector<double> data(num_trials);
	for (r = range[0]; r <= range[1]; r++)
	{
		int cnt_zero = 0;
		for (i = 0; i < num_trials; i++)
		{
			cafe_tree_random_familysize(param->pcafe, r, probability_cache->maxFamilysize);
			data[i] = __cafe_cmd_extinct_count_zero((pTree)param->pcafe);
			cnt_zero += data[i];
		}
		cafe_log(param, "------------------------------------------\n");
		cafe_log(param, "Root size: %d\n", r);
		histogram_set_sparse_data(phist_sim, &data[0], num_trials);
		histogram_merge(phist_accu, phist_sim);
		histogram_print(phist_sim, param->flog);
		if (param->flog != stdout) histogram_print(phist_sim, NULL);
		cafe_log(param, "Sum : %d\n", cnt_zero);
		accu_sum += cnt_zero;
	}

	cafe_log(param, "------------------------------------------\n");
	cafe_log(param, "Total\n", r);
	histogram_print(phist_accu, param->flog);
	if (param->flog != stdout) histogram_print(phist_accu, NULL);
	cafe_log(param, "Sum : %d\n", accu_sum);


	histogram_free(phist_sim);
	histogram_free(phist_accu);
	tree_clear_reg((pTree)param->pcafe);
	return 0;
}

/**
\ingroup Commands
\brief Change the length of an individual branch or all branches, based on a file

*/
int cafe_cmd_branchlength(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

	prereqs(param, REQUIRES_TREE);

	pString pstr = cafe_tree_string_with_id(param->pcafe);
	printf("%s\n", pstr->buf);
	string_free(pstr);
	int err = 0;

	if (tokens.size() == 1)
	{
		cafe_shell_set_branchlength(param, probability_cache->maxFamilysize);
	}
	else if (tokens.size() > 2)
	{
		std::vector<int> lengths(tokens.size() - 1);
		std::transform(tokens.begin() + 1, tokens.end(), lengths.begin(), s_to_i);
		tree_set_branch_lengths(param->pcafe, lengths);
	}
	if (!err)
	{
		cafe_tree_string_print(param->pcafe);
	}
	if (param->pfamily)
	{
		int i;
		for (i = 0; i < param->pfamily->flist->size; i++)
		{
			pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
			pitem->maxlh = -1;
		}
	}
	return err;
}

/**
\ingroup Commands

*/
int cafe_cmd_accuracy(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

	std::string truthfile;
	vector<Argument> args = build_argument_list(tokens);
	for (size_t i = 0; i < args.size(); i++)
	{
		pArgument parg = &args[i];

		if (!strcmp(parg->opt.c_str(), "-i"))
		{
			ostringstream ost;
			copy(parg->argv, parg->argv + parg->argc, std::ostream_iterator<string>(ost, " "));
			truthfile = ost.str();
		}
	}

	int i, j;
	pCafeTree pcafe = param->pcafe;
	double SSE = 0;
	double MSE = 0;
	int nodecnt = 0;
	int errnodecnt = 0;
	if (!truthfile.empty())
	{
		// read in truth data
    ifstream ifst(truthfile);
		pCafeFamily truthfamily = load_gene_families(ifst, '\t', -1);
		if (truthfamily == NULL) {
			fprintf(stderr, "failed to read in true values %s\n", truthfile.c_str());
			return -1;
		}
		pCafeTree truthtree = cafe_tree_copy(pcafe);
		// set parameters
		if (truthtree)
		{
			cafe_family_set_species_index(truthfamily, truthtree);
			// compare inferred vs. truth
			for (i = 0; i< param->pfamily->flist->size; i++)
			{
        pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
        cafe_family_set_truesize_with_family(truthfamily, i, truthtree);
				cafe_family_set_size(param->pfamily, pitem, pcafe);
				if (param->posterior) {
					cafe_tree_viterbi_posterior(pcafe, param);
				}
				else {
					cafe_tree_viterbi(pcafe);
				}
				// inner node SSE
				for (j = 1; j<pcafe->super.nlist->size; j = j + 2) {
					int error = ((pCafeNode)truthtree->super.nlist->array[j])->familysize - ((pCafeNode)pcafe->super.nlist->array[j])->familysize;
					SSE += pow(error, 2);
					if (error != 0) { errnodecnt++; }
					nodecnt++;
				}
			}
			MSE = SSE / nodecnt;
		}

	}
	cafe_log(param, "ancestral reconstruction SSE %f MSE %f totalnodes %d errornodes %d\n", SSE, MSE, nodecnt, errnodecnt);
	return 0;
}

/**
\ingroup Commands
\brief Write gains and losses
*
*/
int cafe_cmd_gainloss(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;
	prereqs(param, REQUIRES_FAMILY | REQUIRES_TREE | REQUIRES_LAMBDA);

	if (globals.viterbi->viterbiNodeFamilysizes.empty())
	{
		if (ConditionalDistribution::matrix.empty())
		{
            cafe_shell_set_lambdas(param, param->input.parameters);
            reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, &param->family_size);
			ConditionalDistribution::reset(param->pcafe, &param->family_size, param->num_threads, globals.num_random_samples);
		}
		pArrayList cd = ConditionalDistribution::to_arraylist();
		cafe_viterbi(globals, *globals.viterbi, &ConditionalDistribution::matrix);
		arraylist_free(cd, NULL);
	}

	string name = tokens[1] + ".gs";
	ofstream ofst(name.c_str());
	//FILE* fout = fopen(name.c_str(), "w");

	pCafeTree pcafe = param->pcafe;
	pCafeTree psum = cafe_tree_copy(pcafe);

	clear_tree_viterbis(psum);

	int totalsum = 0;
	int fsize = param->pfamily->flist->size;

	for (int i = 0; i < fsize; i++)
	{
    pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
    cafe_family_set_size(param->pfamily, pitem, pcafe);
    globals.viterbi->set_node_familysize(pcafe, pitem);
		totalsum += write_family_gainloss(ofst, pitem->id, pcafe, psum);
	}
	ofst << "SUM\t" << totalsum << "\t";
	pString pstr = phylogeny_string((pTree)psum, __cafe_tree_string_sum_gainloss);
	ofst << pstr->buf << "\n";
	string_free(pstr);

	cafe_tree_free(psum);

	return 0;
}

/**
\ingroup Commands
\brief Logs various pieces of information about the application state
*
*/
int cafe_cmd_print_param(Globals& globals, std::vector<std::string> tokens)
{
	log_buffer buf(&globals.param);
	ostream ost(&buf);
	log_param_values(ost, globals);
	return 0;
}

/**
\ingroup Commands
\brief Cross validation by family

Arguments: -fold
*/
int cafe_cmd_cvfamily(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

	prereqs(param, REQUIRES_FAMILY | REQUIRES_TREE | REQUIRES_LAMBDA);

	int cv_fold = 0;
	vector<Argument> args = build_argument_list(tokens);
	for (size_t i = 0; i < args.size(); i++)
	{
		pArgument parg = &args[i];

		if (!strcmp(parg->opt.c_str(), "-fold"))
		{
			sscanf(parg->argv[0], "%d", &cv_fold);
		}
	}

	double MSE_allfolds = 0;
	pCafeFamily pcafe_original = param->pfamily;

	if (tokens.size() < 2)
	{
		throw std::runtime_error("Usage(cvfamily): cvfamily -fold <num>\n");
	}

	// set up the training-validation set
	cafe_family_split_cvfiles_byfamily(param, cv_fold);

	log_buffer buf(&globals.param);
	ostream log(&buf);
	//
	for (int i = 0; i<cv_fold; i++)
	{
		char trainfile[STRING_BUF_SIZE];
		char queryfile[STRING_BUF_SIZE];
		char validatefile[STRING_BUF_SIZE];
		sprintf(trainfile, "%s.%d.train", param->str_fdata->buf, i + 1);
		sprintf(queryfile, "%s.%d.query", param->str_fdata->buf, i + 1);
		sprintf(validatefile, "%s.%d.valid", param->str_fdata->buf, i + 1);

		// read in training data
    ifstream ifst(trainfile);
		pCafeFamily tmpfamily = load_gene_families(ifst, '\t', -1);
		if (tmpfamily == NULL) {
			fprintf(stderr, "failed to read in training data %s\n", trainfile);
			return -1;
		}
		param->pfamily = tmpfamily;

		set_range_from_family(&param->family_size, param->pfamily);
		if (param->pcafe)
		{
			cafe_tree_set_parameters(param->pcafe, &param->family_size, 0);
			cafe_family_set_species_index(param->pfamily, param->pcafe);
		}
		// re-train 
		if (param->num_mus > 0) {
			best_lambda_mu_by_fminsearch(param, param->num_lambdas, param->num_mus, param->parameterized_k_value, log);
		}
		else {
			cafe_best_lambda_by_fminsearch(param, param->num_lambdas, param->parameterized_k_value);
		}

		//cross-validate
		double MSE = globals.validator->validate_by_family(&globals.param, queryfile, validatefile, "MSE");
		MSE_allfolds += MSE;
		cafe_log(param, "MSE fold %d %f\n", i + 1, MSE);

		cafe_family_free(tmpfamily);
	}
	MSE_allfolds = MSE_allfolds / cv_fold;
	cafe_log(param, "MSE all folds %f\n", MSE_allfolds);

	//re-load the original family file
	param->pfamily = pcafe_original;
	set_range_from_family(&param->family_size, param->pfamily);

	if (param->pcafe)
	{
		cafe_tree_set_parameters(param->pcafe, &param->family_size, 0);
		cafe_family_set_species_index(param->pfamily, param->pcafe);
	}
	// re-train 
	if (param->num_mus > 0) {
		best_lambda_mu_by_fminsearch(param, param->num_lambdas, param->num_mus, param->parameterized_k_value, log);
	}
	else {
		cafe_best_lambda_by_fminsearch(param, param->num_lambdas, param->parameterized_k_value);
	}

	// remove training-validation set
  globals.validator->clean_by_family(param->str_fdata->buf, cv_fold);
	return 0;
}

/**
\ingroup Commands
\brief Cross validation by species

*/
int cafe_cmd_cvspecies(Globals& globals, std::vector<std::string> tokens)
{
	log_buffer buf(&globals.param);
	ostream log(&buf);
	pCafeParam param = &globals.param;

	int i;
	prereqs(param, REQUIRES_FAMILY | REQUIRES_TREE | REQUIRES_LAMBDA);

	double MSE_allspecies = 0;
	pCafeFamily pcafe_original = param->pfamily;
	int num_species_original = param->pfamily->num_species;
	char** species_names_original = param->pfamily->species;

	if (tokens.size() < 2)
	{
		// set up the training-validation set
		cafe_family_split_cvfiles_byspecies(param);

		for (i = 0; i<num_species_original; i++) {
			char trainfile[STRING_BUF_SIZE];
			char validatefile[STRING_BUF_SIZE];
			sprintf(trainfile, "%s.%s.train", param->str_fdata->buf, species_names_original[i]);
			sprintf(validatefile, "%s.%s.valid", param->str_fdata->buf, species_names_original[i]);

			// read in training data
      std::ifstream ifst(trainfile);
			pCafeFamily tmpfamily = load_gene_families(ifst, '\t', -1);
			if (tmpfamily == NULL) {
				fprintf(stderr, "failed to read in training data %s\n", trainfile);
				fprintf(stderr, "did you load the family data with the cross-validation option (load -i <familyfile> -cv)?\n");
				return -1;
			}
			param->pfamily = tmpfamily;

			set_range_from_family(&param->family_size, param->pfamily);

			if (param->pcafe)
			{
				cafe_tree_set_parameters(param->pcafe, &param->family_size, 0);
				cafe_family_set_species_index(param->pfamily, param->pcafe);
			}
			// re-train 
			if (param->num_mus > 0) {
				best_lambda_mu_by_fminsearch(param, param->num_lambdas, param->num_mus, param->parameterized_k_value, log);
			}
			else {
				cafe_best_lambda_by_fminsearch(param, param->num_lambdas, param->parameterized_k_value);
			}

			//cross-validate
			double MSE = globals.validator->validate_by_species(param, validatefile, "MSE");
			MSE_allspecies += MSE;
			cafe_log(param, "MSE %s %f\n", globals.validator->get_species_name().c_str(), MSE);

			cafe_family_free(tmpfamily);
		}
		MSE_allspecies = MSE_allspecies / num_species_original;
		cafe_log(param, "MSE all species %f\n", MSE_allspecies);

		//re-load the original family file
		param->pfamily = pcafe_original;
		set_range_from_family(&param->family_size, param->pfamily);

		if (param->pcafe)
		{
			cafe_tree_set_parameters(param->pcafe, &param->family_size, 0);
			cafe_family_set_species_index(param->pfamily, param->pcafe);
		}
		// re-train 
		if (param->num_mus > 0) {
			best_lambda_mu_by_fminsearch(param, param->num_lambdas, param->num_mus, param->parameterized_k_value, log);
		}
		else {
			cafe_best_lambda_by_fminsearch(param, param->num_lambdas, param->parameterized_k_value);
		}

		// remove training-validation set
    globals.validator->clean_by_species(param->str_fdata->buf);
	}
	else
	{
		string validate_file = get_input_file(tokens);
    globals.validator->validate_by_species(param, validate_file.c_str(), "MSE");
	}
	return 0;
}

/**
\ingroup Commands
\brief Runs a simulation to determine how many families are likely to have gone extinct

Arguments: –t parameter gives the number of trials. Logs the total number of extinct families?
*/
int cafe_cmd_extinct(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

	prereqs(param, REQUIRES_FAMILY | REQUIRES_LAMBDA | REQUIRES_TREE);

	int familysize = param->pfamily->flist->size;
	pCafeTree pcafe = param->pcafe;

	int num_trials = 1000;

	cafe_log(param, ">> Data and Simulation for extinction:\n");

	vector<Argument> args = build_argument_list(tokens);
	if (args.size() > 0 && !strcmp(args[0].opt.c_str(), "-t"))
	{
		sscanf(args[0].argv[0], "%d", &num_trials);
	}

	cafe_log(param, "# trials: %d\n", num_trials);

	int i, j;
	int total_extinct = 0;

	fprintf(stderr, "Data ...\n");

	roots roots;
	run_viterbi_sim(pcafe, param->pfamily, roots);

	int maxsize = init_histograms(pcafe->rfsize, roots, familysize);

	pHistogram phist_tmp = histogram_new(NULL, 0, 0);

	double* data = (double*)memory_new(maxsize, sizeof(double));

	cafe_log(&globals.param, "*******************  DATA  *************************\n");

	for (i = 1; i <= pcafe->rfsize; i++)
	{
		if (roots.num[i])
		{
			int k = 0;
			int sum = 0;
			for (j = 0; j < familysize; j++)
			{
				if (roots.size[j] == i)
				{
					data[k++] = roots.extinct[j];
					sum += roots.extinct[j];
				}
			}
			histogram_set_sparse_data(roots.phist_data[i], data, k);
			cafe_log(param, "--------------------------------\n");
			cafe_log(param, "Root Size: %d\n", i);
			ostringstream ost;
			ost << *roots.phist_data[i];
			cafe_log(param, ost.str().c_str());

			//histogram_print(roots.phist_data[i], );
			if (param->flog != stdout) histogram_print(roots.phist_data[i], NULL);
			cafe_log(param, "Extinct: %d\n", sum);
		}
	}

	fprintf(stderr, "Begin simulation...\n");
	cafe_log(param, "******************* SIMULATION **********************\n");

	pHistogram** phist_sim_n = (pHistogram**)memory_new_2dim(num_trials, pcafe->rfsize + 1, sizeof(pHistogram));

	double* simdata = (double*)memory_new(familysize, sizeof(double));
	int t;
	for (t = 0; t < num_trials; t++)
	{
		if (t % 100 == 0 && t != 0)
		{
			fprintf(stderr, "\t%d...\n", t);
		}
		int f = 0;
		for (i = 1; i <= pcafe->rfsize; i++)
		{
			if (roots.num[i])
			{
				int k = 0;
				int sum = 0;
				for (j = 0; j < roots.num[i]; j++)
				{
					cafe_tree_random_familysize(pcafe, i, probability_cache->maxFamilysize);
					simdata[f] = __cafe_cmd_extinct_count_zero((pTree)pcafe);
					data[k++] = simdata[f];
					sum += simdata[f];
					f++;
				}
				roots.avg_extinct[i] += sum;
				histogram_set_sparse_data(phist_tmp, data, k);
				histogram_merge(roots.phist_sim[i], phist_tmp);
				phist_sim_n[t][i] = histogram_new(NULL, 0, 0);
				histogram_merge(phist_sim_n[t][i], phist_tmp);
			}
		}
		histogram_set_sparse_data(phist_tmp, simdata, familysize);
		histogram_merge(roots.phist_sim[0], phist_tmp);
		phist_sim_n[t][0] = histogram_new(NULL, 0, 0);
		histogram_merge(phist_sim_n[t][0], phist_tmp);
	}

	double* cnt = (double*)memory_new(num_trials, sizeof(double));
	double avg_total_extinct = 0;
	for (i = 1; i <= pcafe->rfsize; i++)
	{
		if (roots.phist_sim[i])
		{
			cafe_log(param, "--------------------------------\n");
			cafe_log(param, "Root Size: %d\n", i);
			__hg_print_sim_extinct(phist_sim_n, roots.phist_sim, i, phist_tmp, cnt, num_trials);
			roots.avg_extinct[i] /= num_trials;
			avg_total_extinct += roots.avg_extinct[i];
		}
	}



	cafe_log(param, "*******************  ALL *************************\n");
	cafe_log(param, ">> DATA\n");
	histogram_print(roots.phist_data[0], param->flog);
	if (param->flog != stdout) histogram_print(roots.phist_data[0], NULL);
	cafe_log(param, "Total Extinct: %d\n", total_extinct);
	cafe_log(param, ">> SIMULATION\n");
	__hg_print_sim_extinct(phist_sim_n, roots.phist_sim, 0, phist_tmp, cnt, num_trials);

	memory_free(cnt);
	cnt = NULL;

	for (i = 0; i < pcafe->rfsize; i++)
	{
		if (roots.phist_data[i])
		{
			histogram_free(roots.phist_data[i]);
			histogram_free(roots.phist_sim[i]);
			for (t = 0; t < num_trials; t++)
			{
				histogram_free(phist_sim_n[t][i]);
			}
		}
	}

	histogram_free(phist_tmp);
	memory_free(roots.phist_data);
	roots.phist_data = NULL;
	memory_free(roots.phist_sim);
	roots.phist_sim = NULL;
	memory_free_2dim((void**)phist_sim_n, num_trials, 0, NULL);

	memory_free(data);
	data = NULL;

	memory_free(simdata);
	simdata = NULL;
	return 0;
}

/**
\ingroup Commands
\brief Allows modifications of family data
*
* Takes four arguments: --idx, -id, -add, and -filter
*/
int cafe_cmd_family(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

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

		log_buffer buf(&globals.param);
		ostream ost(&buf);

		log_param_values(ost, globals);
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
		if (param->pcafe && probability_cache) 
            viterbi_family_print(param->pcafe, param->pfamily, idx);
	}

	return pitem ? 0 : -1;
}

/**
\ingroup Commands
\brief Load program state from a file
*
*/
int cafe_cmd_retrieve(Globals& globals, std::vector<std::string> tokens)
{
	globals.Clear(1);
	if (cafe_report_retrieve_data(tokens[1].c_str(), &globals.param, *globals.viterbi) == -1)
	{
		return -1;
	}
	return 0;
}

/**
\ingroup Commands
\brief Save
*
*/
int cafe_cmd_save(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

	if (tokens.size() != 2)
	{
		throw std::runtime_error("Usage(save): save filename");
	}
	ofstream ofst(tokens[1].c_str());
	if (!ofst)
	{
		throw io_error("save", tokens[1], true);
	}
	write_family(ofst, param->pfamily);

	return 0;
}

void viterbi_write(ostream& ost, pCafeTree pcafe, pCafeFamily pfamily)
{
	pCafeFamilyItem pitem;
	double score = 0;
	for (int i = 0; i < pfamily->flist->size; i++)
	{
		pitem = (pCafeFamilyItem)pfamily->flist->array[i];
		pitem->maxlh = -1;
		cafe_family_set_size_with_family(pfamily, i, pcafe);
		compute_tree_likelihoods(pcafe);
		int ridx = __maxidx(((pCafeNode)pcafe->super.root)->likelihoods, pcafe->rfsize) + pcafe->rootfamilysizes[0];
		double mlh = __max(((pCafeNode)pcafe->super.root)->likelihoods, pcafe->rfsize);
		score += log(mlh);
		cafe_tree_viterbi(pcafe);
		pString pstr = cafe_tree_string(pcafe);
		ost << pitem->id << "\t" << mlh << "\t" << pstr->buf << "\t" << ridx << "\n";
		string_free(pstr);
	}
	ost << "Score: " << score << "\n";
}


/**
\ingroup Commands
\brief Viterbi
*
*/
int cafe_cmd_viterbi(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;

	prereqs(param, REQUIRES_TREE | REQUIRES_LAMBDA);
	struct viterbi_args args = get_viterbi_arguments(build_argument_list(tokens));


	if (args.all)
	{
		cafe_tree_set_parameters(param->pcafe, &param->family_size, 0);
		prereqs(param, REQUIRES_FAMILY);
		ostream* fp = &cout;
		ofstream fout;
		if (!args.file.empty())
		{
			fout.open(args.file.c_str());
			if (fout.fail())
			{
				throw io_error("viterbi", args.file, true);
			}
			fp = &fout;
		}
		viterbi_write(*fp, param->pcafe, param->pfamily);
	}
	else if (!args.item_id.empty())
	{
		int idx = cafe_family_get_index(param->pfamily, args.item_id.c_str());
		if (idx == -1)
		{
			ostringstream ost;
			ost << "ERROR(viterbi): " << args.item_id << " not found";
			throw std::runtime_error(ost.str());
		}
		viterbi_family_print(param->pcafe, param->pfamily, idx);
	}
	else if (args.idx >= 0)
	{
		if (args.idx > param->pfamily->flist->size)
		{
			ostringstream ost;
			ost << "ERROR(viterbi): Out of range[0~" << param->pfamily->flist->size << "]: " << args.idx;
			throw std::runtime_error(ost.str());
		}
		viterbi_family_print(param->pcafe, param->pfamily, args.idx);
	}
	else if (tokens.size() == 1)
	{
		viterbi_print(param->pcafe, set_family_size_interactive(param->pcafe));
	}
	else
	{
		tokens.erase(tokens.begin());
		viterbi_print(param->pcafe, cafe_shell_parse_familysize((pTree)param->pcafe, tokens));
	}

	return 0;
}


#endif
