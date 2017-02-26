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

extern "C" {
	extern pCafeParam cafe_param;
	void cafe_shell_set_lambda(pCafeParam param, double* parameters);
	int __cafe_cmd_lambda_tree(pArgument parg);
}

using namespace std;

int get_doubles_array(double** loc, pArgument parg)
{
	if (*loc) memory_free(*loc);
	*loc = (double*)memory_new(parg->argc, sizeof(double));
	double *dlist = *loc;
	for (int j = 0; j < parg->argc; j++)
	{
		sscanf(parg->argv[j], "%lf", &dlist[j]);
	}
	return parg->argc;
}

bool is_out(Argument arg)
{
	return !strcmp(arg.opt, "-o");
}

lambda_args get_arguments(vector<Argument> pargs)
{
	lambda_args result;

	for (size_t i = 0; i < pargs.size(); i++)
	{
		pArgument parg = &pargs[i];

		// Search for whole family 
		if (!strcmp(parg->opt, "-s"))
		{
			result.search = true;
		}
		else if (!strcmp(parg->opt, "-checkconv"))
		{
			result.tmp_param.checkconv = 1;
		}
		else if (!strcmp(parg->opt, "-t"))
		{
			result.bdone = __cafe_cmd_lambda_tree(parg);
			if (result.bdone < 0) {
				throw std::exception();
			}
			pString pstr = phylogeny_string(cafe_param->lambda_tree, NULL);
			cafe_log(cafe_param, "Lambda Tree: %s\n", pstr->buf);
			string_free(pstr);
			result.tmp_param.lambda_tree = cafe_param->lambda_tree;
			result.tmp_param.num_lambdas = cafe_param->num_lambdas;

		}
		else if (!strcmp(parg->opt, "-v"))
		{
			sscanf(parg->argv[0], "%lf", &result.vlambda);
			result.lambda_type = SINGLE_LAMBDA;
		}
		else if (!strcmp(parg->opt, "-l"))
		{
			result.tmp_param.num_params += get_doubles_array(&result.tmp_param.lambda, parg);
			result.lambda_type = MULTIPLE_LAMBDAS;

		}
		else if (!strcmp(parg->opt, "-p"))
		{
			// TODO: Are the l and p parameters really supposed to be added together into num_params?
			result.tmp_param.num_params += get_doubles_array(&result.tmp_param.k_weights, parg);
		}
		else if (!strcmp(parg->opt, "-r"))
		{
			result.dist = *parg;
			vector<Argument>::iterator it = find_if(pargs.begin(), pargs.end(), is_out);
			if (it != pargs.end())
				result.out = *it;
		}
		else if (!strcmp(parg->opt, "-e"))
		{
			vector<Argument>::iterator it = find_if(pargs.begin(), pargs.end(), is_out);
			if (it != pargs.end())
				result.out = *it;
			result.write_files = true;
			result.each = true;
		}
		else if (!strcmp(parg->opt, "-k"))
		{
			sscanf(parg->argv[0], "%d", &result.tmp_param.parameterized_k_value);
		}
		else if (!strcmp(parg->opt, "-f"))
		{
			result.tmp_param.fixcluster0 = 1;
		}
	}

	return result;
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

void write_lambda_distribution(pArgument parg, FILE* fp)
{
	double** range = (double**)memory_new_2dim(parg->argc, 3, sizeof(double));
	int j;
	for (j = 0; j < parg->argc; j++)
	{
		sscanf(parg->argv[j], "%lf:%lf:%lf", &range[j][0], &range[j][1], &range[j][2]);
		cafe_log(cafe_param, "%dst Distribution: %lf : %lf : %lf\n", j + 1, range[j][0], range[j][1], range[j][2]);
	}
	cafe_param->num_lambdas = parg->argc;
	pGMatrix pgm = cafe_lambda_distribution(cafe_param, parg->argc, range);
	if (fp)
	{
		int k;
		int* idx = (int*)memory_new(parg->argc, sizeof(int));
		for (j = 0; j < pgm->num_elements; j++)
		{
			gmatrix_dim_index(pgm, j, idx);
			fprintf(fp, "%lf", idx[0] * range[0][1] + range[0][0]);
			for (k = 1; k < parg->argc; k++)
			{
				fprintf(fp, "\t%lf", idx[k] * range[k][1] + range[k][0]);
			}
			fprintf(fp, "\t%lf\n", gmatrix_double_get_with_index(pgm, j));
		}
		fclose(fp);
		memory_free(idx);
		idx = NULL;
	}
	memory_free_2dim((void**)range, parg->argc, 0, NULL);
	gmatrix_free(pgm);
}


void validate_lambda_count(int expected, int actual, pTree pTree, int k_value)
{
	// check if the numbers of lambdas and proportions put in matches the number of parameters
	if (expected != actual) {
		ostringstream ost;
		ost << "ERROR(lambda): Number of parameters not correct. \n";
		ost << "The number of -l lambdas are " << actual << "they need to be " << expected << "\n";
		if (pTree || k_value > 0)
			ost << "Based on";
		if (pTree)
		{
			pString pstr = phylogeny_string(pTree, NULL);
			ost << " the tree " << pstr->buf;
			string_free(pstr);
		}
		if (k_value > 0)
		{
			if (pTree) ost << " and ";
			ost << "the -k clusters " << k_value << ".\n";
		}
		throw std::runtime_error(ost.str().c_str());
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
	
	int i,j;

	if(!param->pcafe)
	{
		throw std::runtime_error("ERROR(lambda): You did not specify tree: command 'tree'\n" );
	}

	pCafeTree pcafe = param->pcafe;
	vector<Argument> pargs = build_argument_list(tokens);

	globals.Prepare();

	lambda_args params = get_arguments(pargs);

	if (params.lambda_type == SINGLE_LAMBDA && params.vlambda > 0 )
	{
		set_all_lambdas(param, params.vlambda);
	}
	if ( params.dist.opt )
	{
		FILE* fp = NULL;
		if( params.out.opt && (fp = fopen(params.out.argv[0],"w") ) == NULL )
		{
			fprintf( stderr, "ERROR(lambda): Cannot open file: %s\n", params.out.argv[0] );
			return -1;
		}
		param->posterior = params.tmp_param.posterior;
		if (param->posterior) {
			// set rootsize prior based on leaf size
			cafe_set_prior_rfsize_empirical(param);
		}		
		param->num_params = param->num_lambdas;
        
		if( param->parameters ) memory_free(param->parameters);
		param->parameters = NULL;
		param->parameters = (double*)memory_new(param->num_params, sizeof(double));

		write_lambda_distribution(&params.dist, fp);
		params.bdone = 1;
		fclose(fp);
	}

	if (params.bdone )
	{
		if (params.bdone ) return 0;
	}

	// copy parameters collected to param based on the combination of options.
	{
		param->posterior = params.tmp_param.posterior;
		if (param->posterior) {
			// set rootsize prior based on leaf size
			cafe_set_prior_rfsize_empirical(param);
		}		
		// search or set
		if (params.search) {
            // prepare parameters
			if (params.tmp_param.lambda_tree != NULL) {
				// param->num_lambdas determined by lambda tree.
				if (params.tmp_param.parameterized_k_value > 0) {
					param->parameterized_k_value = params.tmp_param.parameterized_k_value;
					param->fixcluster0 = params.tmp_param.fixcluster0;
					param->num_params = (params.tmp_param.num_lambdas*(params.tmp_param.parameterized_k_value- params.tmp_param.fixcluster0))+(params.tmp_param.parameterized_k_value-1);

					if( param->parameters ) memory_free(param->parameters);
					param->parameters = NULL;
					param->parameters = (double*)memory_new(param->num_params, sizeof(double));
					if (param->k_weights) { memory_free(param->k_weights);}
					param->k_weights = NULL;
					param->k_weights = (double*) memory_new(param->parameterized_k_value, sizeof(double));

					//param->lambda = param->parameters;
				}
				else {	// search whole dataset branch specific
					param->num_params = param->num_lambdas;
					
					if( param->parameters ) memory_free(param->parameters);
					param->parameters = NULL;
					param->parameters = (double*)memory_new(param->num_params, sizeof(double));
					
					//param->lambda = param->parameters;
				}
			}
			else {
				param->num_lambdas = params.tmp_param.num_lambdas = 1;
				if (params.tmp_param.parameterized_k_value > 0) {
					param->parameterized_k_value = params.tmp_param.parameterized_k_value;
					param->fixcluster0 = params.tmp_param.fixcluster0;
					param->num_params = (params.tmp_param.num_lambdas*(params.tmp_param.parameterized_k_value- params.tmp_param.fixcluster0))+(params.tmp_param.parameterized_k_value-1);
					
					if( param->parameters ) memory_free(param->parameters);
					param->parameters = NULL;
					param->parameters = (double*)memory_new(param->num_params, sizeof(double));
					if (param->k_weights) { memory_free(param->k_weights);}
					param->k_weights = NULL;
					param->k_weights = (double*) memory_new(param->parameterized_k_value, sizeof(double));
					
					//param->lambda = param->parameters;
				}
				else {	// search whole dataset whole tree
					param->num_params = param->num_lambdas;
					
					if( param->parameters ) memory_free(param->parameters);
					param->parameters = NULL;
					param->parameters = (double*)memory_new(param->num_params, sizeof(double));
					
					//param->lambda = param->parameters;
				}
			}
			// search
			if (params.tmp_param.checkconv) { param->checkconv = 1; }
			if (params.each )
			{
				cafe_each_best_lambda_by_fminsearch(param,param->num_lambdas);
			}
			else
			{
				cafe_best_lambda_by_fminsearch(param, param->num_lambdas, param->parameterized_k_value);
			}
		}
		else {
			if (params.tmp_param.lambda_tree != NULL) {
				// param->num_lambdas determined by lambda tree.
				if (params.tmp_param.parameterized_k_value > 0) {	// search clustered branch specific
					param->parameterized_k_value = params.tmp_param.parameterized_k_value;
					param->fixcluster0 = params.tmp_param.fixcluster0;
					param->num_params = (params.tmp_param.num_lambdas*(params.tmp_param.parameterized_k_value- params.tmp_param.fixcluster0))+(params.tmp_param.parameterized_k_value-1);

					validate_lambda_count(param->num_params, params.tmp_param.num_params, params.tmp_param.lambda_tree, param->parameterized_k_value);
					
					// copy user input into parameters
					if( param->parameters ) memory_free(param->parameters);
					param->parameters = NULL;
					param->parameters = (double*)memory_new(param->num_params, sizeof(double));
					memcpy(param->parameters,params.tmp_param.lambda, sizeof(double)*params.tmp_param.num_lambdas*(params.tmp_param.parameterized_k_value-params.tmp_param.fixcluster0));
					memcpy(&param->parameters[param->num_lambdas*(param->parameterized_k_value-params.tmp_param.fixcluster0)], params.tmp_param.k_weights, sizeof(double)*(params.tmp_param.parameterized_k_value-1));
					// prepare space for k_weights
					if ( param->k_weights ) memory_free(param->k_weights);
					param->k_weights = NULL;
					param->k_weights = (double*)memory_new(param->parameterized_k_value-1, sizeof(double) );										
					
				}
				else {	// search whole dataset branch specific
					param->num_params = param->num_lambdas;
					
					validate_lambda_count(param->num_params, params.tmp_param.num_params, params.tmp_param.lambda_tree, -1);
					
					// copy user input into parameters
					if( param->parameters ) memory_free(param->parameters);
					param->parameters = NULL;
					param->parameters = (double*)memory_new(param->num_params, sizeof(double));
					memcpy(param->parameters,params.tmp_param.lambda, sizeof(double)*params.tmp_param.num_lambdas);
				}
			}
			else {
				param->num_lambdas = params.tmp_param.num_lambdas = 1;	
				if (params.tmp_param.parameterized_k_value > 0) {				// search clustered whole tree
					param->parameterized_k_value = params.tmp_param.parameterized_k_value;
					param->fixcluster0 = params.tmp_param.fixcluster0;
					param->num_params = (param->num_lambdas*(param->parameterized_k_value-param->fixcluster0))+(param->parameterized_k_value-1);
					
					validate_lambda_count(param->num_params, params.tmp_param.num_params, NULL, param->parameterized_k_value);
					
					// copy user input into parameters
					if( param->parameters ) memory_free(param->parameters);
					param->parameters = NULL;
					param->parameters = (double*)memory_new(param->num_params, sizeof(double));
					memcpy(param->parameters,params.tmp_param.lambda, sizeof(double)*params.tmp_param.num_lambdas*(params.tmp_param.parameterized_k_value-params.tmp_param.fixcluster0));
					memcpy(&param->parameters[param->num_lambdas*(param->parameterized_k_value-params.tmp_param.fixcluster0)], params.tmp_param.k_weights, sizeof(double)*(params.tmp_param.parameterized_k_value-1));
					// prepare space for k_weights
					if ( param->k_weights ) memory_free(param->k_weights);
					param->k_weights = NULL;
					param->k_weights = (double*)memory_new(param->parameterized_k_value-1, sizeof(double) );										
					
				}
				else {	// search whole dataset whole tree
					param->num_params = param->num_lambdas;
					
					validate_lambda_count(param->num_params, params.tmp_param.num_params, NULL, -1);

					// copy user input into parameters
					if( param->parameters ) memory_free(param->parameters);
					param->parameters = NULL;
					param->parameters = (double*)memory_new(param->num_params, sizeof(double));
					memcpy(param->parameters,params.tmp_param.lambda, sizeof(double)*params.tmp_param.num_lambdas);
				}
			}
			param->param_set_func(param, param->parameters);
		}
	}
		
	FILE* fpout = stdout;
	FILE* fmpout = NULL;
	FILE* fhttp = NULL;
	if (params.write_files)
	{
		string base_name = params.out.argv[0];
		params.name = base_name + ".lambda";
		if ((fpout = fopen(params.name.c_str(), "w")) == NULL)
		{
			throw std::runtime_error((std::string("Cannot open file: ") + params.name).c_str());
		}
		params.name = base_name + ".mp";
		if ((fmpout = fopen(params.name.c_str(), "w")) == NULL)
		{
			fclose(fpout);
			throw std::runtime_error((std::string("Cannot open file: ") + params.name).c_str());
		}
		params.name = base_name + ".html";
		if ((fhttp = fopen(params.name.c_str(), "w")) == NULL)
		{
			fclose(fpout);
			fclose(fmpout);
			throw std::runtime_error((std::string("Cannot open file: ") + params.name).c_str());
		}
		params.name = base_name;
	}

	// now print output
	if( params.each )
	{
		if (fhttp )
		{
			fprintf(fhttp,"<html>\n<body>\n<table border=1>\n");
		}
		for ( i = 0 ; i < param->pfamily->flist->size ; i++ )
		{
			pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
			cafe_family_set_size(param->pfamily, i, pcafe);
			param->param_set_func(param,pitem->lambda);
			pString pstr = cafe_tree_string_with_familysize_lambda(pcafe);
			for ( j = 0 ; j < param->num_lambdas; j++ )
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
			if (fmpout )
			{
				pstr = cafe_tree_metapost( pcafe, i+1, pitem->id, 6, 6 );
				fprintf(fmpout, "%s\n", pstr->buf );
				string_free(pstr);
			}
		}
		if (fpout != stdout ) fclose(fpout);
		if (fmpout ) fclose(fmpout);
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



