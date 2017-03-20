extern "C" {
#include <family.h>
#include "cafe_shell.h"
#include "cafe.h"
}

#include<mathfunc.h>
#include <exception>
#include <string>
#include <stdexcept>
#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>

#include "lambda.h"
#include "cafe_commands.h"
#include "Globals.h"
#include "log_buffer.h"

extern "C" {
	extern pCafeParam cafe_param;
	void cafe_shell_set_lambda(pCafeParam param, double* parameters);
	int __cafe_cmd_lambda_tree(pArgument parg);
	double __cafe_best_lambda_search(double* plambda, void* args);
}

using namespace std;

bool is_out(Argument arg)
{
	return !strcmp(arg.opt, "-o");
}

void lambda_arg_base::load(vector<Argument> pargs)
{
	for (size_t i = 0; i < pargs.size(); i++)
	{
		pArgument parg = &pargs[i];

		// Search for whole family 
		if (!strcmp(parg->opt, "-s"))
		{
			search = true;
		}
		else if (!strcmp(parg->opt, "-checkconv"))
		{
			checkconv = true;
		}
		else if (!strcmp(parg->opt, "-t"))
		{
			bdone = __cafe_cmd_lambda_tree(parg);
			if (bdone < 0) {
				throw std::exception();
			}
			pString pstr = phylogeny_string(cafe_param->lambda_tree, NULL);
			cafe_log(cafe_param, "Lambda Tree: %s\n", pstr->buf);
			string_free(pstr);
			lambda_tree = cafe_param->lambda_tree;
			lambdas.resize(cafe_param->num_lambdas);

		}
		else if (!strcmp(parg->opt, "-l"))
		{
			get_doubles_array(lambdas, parg);
			num_params += lambdas.size();
			lambda_type = MULTIPLE_LAMBDAS;

		}
		else if (!strcmp(parg->opt, "-p"))
		{
			get_doubles_array(k_weights, parg);
			num_params += k_weights.size();
		}
		else if (!strcmp(parg->opt, "-k"))
		{
			int k = 0;
			sscanf(parg->argv[0], "%d", &k);
			k_weights.resize(k);
		}
		else if (!strcmp(parg->opt, "-f"))
		{
			fixcluster0 = 1;
		}
	}
}

void lambda_args::load(std::vector<Argument> pargs)
{
	lambda_arg_base::load(pargs);
	for (size_t i = 0; i < pargs.size(); i++)
	{
		pArgument parg = &pargs[i];
		if (!strcmp(parg->opt, "-v"))
		{
			sscanf(parg->argv[0], "%lf", &vlambda);
			lambda_type = SINGLE_LAMBDA;
		} 
		else if (!strcmp(parg->opt, "-r"))
		{
			range.resize(parg->argc);
			int j;
			for (j = 0; j < parg->argc; j++)
			{
				sscanf(parg->argv[j], "%lf:%lf:%lf", &range[j].start, &range[j].step, &range[j].end);
			}
		}
		else if (!strcmp(parg->opt, "-e"))
		{
			write_files = true;
			each = true;
		}
		else if (!strcmp(parg->opt, "-o"))
		{
			outfile = parg->argv[0];
		}
	}
}

void set_all_lambdas(pCafeParam param, double value)
{
	if (param->lambda) memory_free(param->lambda);
	param->lambda = NULL;
	param->num_lambdas = param->num_lambdas < 1 ? 1 : param->num_lambdas;
	param->lambda = (double*)memory_new(param->num_lambdas, sizeof(double));
	for (int j = 0; j < param->num_lambdas; j++)
	{
		param->lambda[j] = value;
	}

}

pGMatrix cafe_lambda_distribution(pCafeParam param, const vector<lambda_range>& range)
{
	int numrange = range.size();
	int i, j;
	int* size = (int*)memory_new(numrange, sizeof(int));
	double* plambda = (double*)memory_new(numrange, sizeof(double));
	int* idx = (int*)memory_new(numrange, sizeof(int));
	for (i = 0; i < numrange; i++)
	{
		size[i] = 1 + rint((range[i].end - range[i].start) / range[i].step);
	}
	pGMatrix pgm = gmatrix_double_new(numrange, size);

	for (i = 0; i < pgm->num_elements; i++)
	{
		gmatrix_dim_index(pgm, i, idx);
		for (j = 0; j < numrange; j++)
		{
			plambda[j] = range[j].step * idx[j] + range[j].start;
		}
		double v = -__cafe_best_lambda_search(plambda, (void*)param);
		gmatrix_double_set_with_index(pgm, v, i);
		if (-v > 1e300)
		{
			for (j = 0; j < param->pfamily->flist->size; j++)
			{
				pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[j];
				pitem->maxlh = -1;
			}
		}
	}

	memory_free(idx);
	idx = NULL;
	memory_free(size);
	size = NULL;
	memory_free(plambda);
	plambda = NULL;
	return pgm;
}

void write_lambda_distribution(pCafeParam param, const std::vector<lambda_range>& range, FILE* fp)
{
	param->num_lambdas = range.size();
	pGMatrix pgm = cafe_lambda_distribution(param, range);
	if (fp)
	{
		int* idx = (int*)memory_new(range.size(), sizeof(int));
		for (int j = 0; j < pgm->num_elements; j++)
		{
			gmatrix_dim_index(pgm, j, idx);
			fprintf(fp, "%lf", idx[0] * range[0].step + range[0].start);
			for (size_t k = 1; k < range.size(); k++)
			{
				fprintf(fp, "\t%lf", idx[k] * range[k].step + range[k].start);
			}
			fprintf(fp, "\t%lf\n", gmatrix_double_get_with_index(pgm, j));
		}
		fclose(fp);
		memory_free(idx);
		idx = NULL;
	}
	gmatrix_free(pgm);
}

void initialize_params_and_k_weights(pCafeParam param, int what)
{
	if (what & INIT_PARAMS)
	{
		if (param->parameters) 
			memory_free(param->parameters);
		param->parameters = (double*)memory_new(param->num_params, sizeof(double));
	}
	if (what & INIT_KWEIGHTS)
	{
		if (param->k_weights) 
			memory_free(param->k_weights);
		param->k_weights = (double*)memory_new(param->parameterized_k_value, sizeof(double));
	}
}

void set_parameters(pCafeParam param, lambda_args& params)
{
	param->parameterized_k_value = params.k_weights.size();
	param->fixcluster0 = params.fixcluster0;
	param->num_params = (params.lambdas.size()*(params.k_weights.size() - params.fixcluster0)) + (params.k_weights.size() - 1);

	params.validate_parameter_count(param->num_params);

	initialize_params_and_k_weights(param, INIT_PARAMS | INIT_KWEIGHTS);

	copy(params.lambdas.begin(), params.lambdas.end(), param->parameters);

	int num_k_weights = params.k_weights.size() - 1;
	int first_k_weight = (param->num_lambdas*(param->parameterized_k_value - params.fixcluster0));
	copy(params.k_weights.begin(), params.k_weights.begin() + num_k_weights, param->parameters + first_k_weight);

}

void lambda_set(pCafeParam param, lambda_args& params)
{
	if (params.lambda_tree != NULL) {
		// param->num_lambdas determined by lambda tree.
		if (params.k_weights.size() > 0) {	// search clustered branch specific
			set_parameters(param, params);
		}
		else {	// search whole dataset branch specific
			param->num_params = param->num_lambdas;

			params.validate_parameter_count(param->num_params);

			// copy user input into parameters
			initialize_params_and_k_weights(param, INIT_PARAMS);
			copy(params.lambdas.begin(), params.lambdas.end(), param->parameters);
			//memcpy(param->parameters, params.lambda, sizeof(double)*params.num_lambdas);
		}
	}
	else {
		param->num_lambdas = 1;
		if (params.k_weights.size() > 0) {				// search clustered whole tree
			set_parameters(param, params);
		}
		else {	// search whole dataset whole tree
			param->num_params = param->num_lambdas;

			params.validate_parameter_count(param->num_params);

			// copy user input into parameters
			initialize_params_and_k_weights(param, INIT_PARAMS);
			//memcpy(param->parameters, params.lambda, sizeof(double)*params.num_lambdas);
			copy(params.lambdas.begin(), params.lambdas.begin() + param->num_lambdas, param->parameters);
		}
	}
	param->param_set_func(param, param->parameters);

}

void lambda_search(pCafeParam param, lambda_args& params)
{
	if (params.lambda_tree == NULL)
	{
		param->num_lambdas = 1;
		params.lambdas.resize(1);
	}
	// prepare parameters
	// param->num_lambdas determined by lambda tree.
	if (params.k_weights.size() > 0) {
		param->parameterized_k_value = params.k_weights.size();
		param->fixcluster0 = params.fixcluster0;
		param->num_params = (params.lambdas.size()*(params.k_weights.size() - params.fixcluster0)) + (params.k_weights.size() - 1);
	}
	else {	// search whole dataset branch specific
		param->num_params = param->num_lambdas;
	}

	initialize_params_and_k_weights(param, INIT_PARAMS | (params.k_weights.size() > 0 ? INIT_KWEIGHTS : 0));

	// search
	if (params.checkconv) { param->checkconv = 1; }
	if (params.each)
	{
		cafe_each_best_lambda_by_fminsearch(param, param->num_lambdas);
	}
	else
	{
		cafe_best_lambda_by_fminsearch(param, param->num_lambdas, param->parameterized_k_value);
	}

}

/**
* \brief Find lambda values
*
* Arguments include -s and -t. The -s argument starts a search
* for lambda values maximizing the log likelihood of data for
* all values. Subsequent analayses will automatically use the
* results from the lambda search.
* -t takes the same Newick tree structure as in the tree
* command, excluding branch lengths and subsituting integer
* values from 1 to N taxon names.
* etc.
*/
int cafe_cmd_lambda(Globals& globals, vector<string> tokens)
{
	pCafeParam param = &globals.param;
	
	prereqs(param, REQUIRES_FAMILY | REQUIRES_TREE);

	pCafeTree pcafe = param->pcafe;
	vector<Argument> pargs = build_argument_list(tokens);

	globals.Prepare();

	lambda_args params;
	params.load(pargs);

	if (params.lambda_type == SINGLE_LAMBDA && params.vlambda > 0 )
	{
		set_all_lambdas(param, params.vlambda);
	}
	if (!params.range.empty())
	{
		FILE* fp = NULL;
		if(!params.outfile.empty() && (fp = fopen(params.outfile.c_str(),"w") ) == NULL )
		{
			fprintf(stderr, "ERROR(lambda): Cannot open file: %s\n", params.outfile.c_str());
			return -1;
		}
		param->posterior = 1;
		// set rootsize prior based on leaf size
		cafe_set_prior_rfsize_empirical(param);

		param->num_params = param->num_lambdas;
        
		initialize_params_and_k_weights(param, INIT_PARAMS);

		log_buffer buf(&globals.param);
		ostream ost(&buf);
		for (size_t j = 0; j < params.range.size(); j++)
		{
			ost << j + 1 << "st Distribution: " << params.range[j].start << " : " << params.range[j].step << " : " << params.range[j].end << "\n";
		}
		write_lambda_distribution(&globals.param, params.range, fp);
		params.bdone = 1;
		fclose(fp);
	}

	if (params.bdone )
	{
		if (params.bdone ) return 0;
	}

	// copy parameters collected to param based on the combination of options.
	param->posterior = 1;
	// set rootsize prior based on leaf size
	cafe_set_prior_rfsize_empirical(param);

	// search or set
	if (params.search) {
		lambda_search(param, params);
	}
	else {
		lambda_set(param, params);
	}
		
	FILE* fpout = stdout;
	FILE* fhttp = NULL;
	if (params.write_files)
	{
		params.name = params.outfile + ".lambda";
		if ((fpout = fopen(params.name.c_str(), "w")) == NULL)
		{
			throw std::runtime_error((std::string("Cannot open file: ") + params.name).c_str());
		}
		params.name = params.outfile + ".html";
		if ((fhttp = fopen(params.name.c_str(), "w")) == NULL)
		{
			fclose(fpout);
			throw std::runtime_error((std::string("Cannot open file: ") + params.name).c_str());
		}
		params.name = params.outfile;
	}

	// now print output
	if( params.each )
	{
		if (fhttp )
		{
			fprintf(fhttp,"<html>\n<body>\n<table border=1>\n");
		}
		for (int i = 0 ; i < param->pfamily->flist->size ; i++ )
		{
			pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
			cafe_family_set_size(param->pfamily, i, pcafe);
			param->param_set_func(param,pitem->lambda);
			pString pstr = cafe_tree_string_with_familysize_lambda(pcafe);
			for (int j = 0 ; j < param->num_lambdas; j++ )
			{
				double a = pitem->lambda[j] * param->max_branch_length;
				if ( a >= 0.5 || fabs(a-0.5) < 1e-3 )
				{
					fprintf(fpout, "@@ ");
					break;
				}	
			}
			fprintf(fpout, "%s\t%s\n", pitem->id, pstr->buf );
			if (fhttp )
			{
				fprintf(fhttp,"<tr><td rowspan=2><a href=pdf/%s-%d.pdf>%s</a></td><td>%s</td></tr>\n",
					params.name.c_str(), i+1, pitem->id, pitem->desc ? pitem->desc : "NONE" );
				fprintf(fhttp,"<tr><td>%s</td></tr>\n", pstr->buf );
			}
			string_free(pstr);
		}
		if (fpout != stdout ) fclose(fpout);
		if (fhttp )
		{
			fprintf(fhttp,"</table>\n</body>\n</html>\n");
			fclose(fhttp );
		}
	}
	else
	{
		if ( param->pfamily )
		{
			reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, &param->family_size);
		}
	}
	
	cafe_log(param,"DONE: Lambda Search or setting, for command:\n");

	std::ostringstream cmd;
	std::copy(tokens.begin(), tokens.end(),
		std::ostream_iterator<std::string>(cmd, " "));
	cafe_log(param,"%s\n", cmd.str().c_str());
	
	if (params.search && (param->parameterized_k_value > 0)) {
		// print the cluster memberships
		cafe_family_print_cluster_membership(param);
	}

	return 0;
}



