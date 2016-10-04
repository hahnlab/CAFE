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

extern "C" {
	extern pCafeParam cafe_param;
	void cafe_shell_set_lambda(pCafeParam param, double* parameters);
	int __cafe_cmd_lambda_tree(pArgument parg);
	void __cafe_cmd_lambda_distribution(pArgument parg, FILE* fp);
}

using namespace std;

vector<Argument> lambda_build_argument(vector<string> tokens)
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
			result.blambda = 2;
		}
		else if (!strcmp(parg->opt, "-l"))
		{
			result.tmp_param.num_params += get_doubles_array(&result.tmp_param.lambda, parg);
			result.blambda = 1;

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
			sscanf(parg->argv[0], "%d", &result.tmp_param.k);
		}
		else if (!strcmp(parg->opt, "-f"))
		{
			result.tmp_param.fixcluster0 = 1;
		}
	}

	return result;
}

void prepare_cafe_param(pCafeParam param)
{
	param->lambda = NULL;
	param->mu = NULL;
	if (param->lambda_tree) {
		phylogeny_free(param->lambda_tree);
		param->lambda_tree = NULL;
	}
	if (param->mu_tree)
	{
		phylogeny_free(param->mu_tree);
		param->mu_tree = NULL;
	}
	param->num_lambdas = -1;
	param->num_mus = -1;
	param->k = 0;
	param->param_set_func = cafe_shell_set_lambda;

	if (param->pfamily == NULL || param->pcafe == NULL)
	{
		throw runtime_error("ERROR(lambda): Please load family (\"load\") and cafe tree (\"tree\") before running \"lambda\" command.");
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
int cafe_cmd_lambda(vector<string> tokens)
{
	try {
		int i,j;

		if(!cafe_param->pcafe)
		{
			fprintf( stderr, "ERROR(lambda): You did not specify tree: command 'tree'\n" );
			return -1;
		}
		pCafeTree pcafe = cafe_param->pcafe;
		vector<Argument> pargs = lambda_build_argument(tokens);

		prepare_cafe_param(cafe_param);

		lambda_args params = get_arguments(pargs);

		// now set or search
		if (params.blambda == 2 && params.vlambda > 0 )
		{
			if ( cafe_param->lambda ) memory_free(cafe_param->lambda);
			cafe_param->lambda = NULL;
			cafe_param->num_lambdas = cafe_param->num_lambdas < 1 ? 1 : cafe_param->num_lambdas;
			cafe_param->lambda = (double*)memory_new( cafe_param->num_lambdas, sizeof(double) );
			for( j = 0 ; j < cafe_param->num_lambdas; j++ )
			{
				cafe_param->lambda[j] = params.vlambda;
			}
		}
		if ( params.dist.opt )
		{
			FILE* fp = NULL;
			if( params.out.opt && (fp = fopen(params.out.argv[0],"w") ) == NULL )
			{
				fprintf( stderr, "ERROR(lambda): Cannot open file: %s\n", params.out.argv[0] );
				return -1;
			}
			cafe_param->posterior = params.tmp_param.posterior;
			if (cafe_param->posterior) {
				// set rootsize prior based on leaf size
				cafe_set_prior_rfsize_empirical(cafe_param);
			}		
			cafe_param->num_params = cafe_param->num_lambdas;
        
			if( cafe_param->parameters ) memory_free(cafe_param->parameters);
			cafe_param->parameters = NULL;
			cafe_param->parameters = (double*)memory_new(cafe_param->num_params, sizeof(double));

			__cafe_cmd_lambda_distribution(&params.dist, fp);
			params.bdone = 1;
			fclose(fp);
		}

	if (params.bdone )
	{
		if (params.bdone ) return 0;
	}

	// copy parameters collected to cafe_param based on the combination of options.
	{
		cafe_param->posterior = params.tmp_param.posterior;
		if (cafe_param->posterior) {
			// set rootsize prior based on leaf size
			cafe_set_prior_rfsize_empirical(cafe_param);
		}		
		// search or set
		if (params.search) {
            // prepare parameters
			if (params.tmp_param.lambda_tree != NULL) {
				// cafe_param->num_lambdas determined by lambda tree.
				if (params.tmp_param.k > 0) {
					cafe_param->k = params.tmp_param.k;
					cafe_param->fixcluster0 = params.tmp_param.fixcluster0;
					cafe_param->num_params = (params.tmp_param.num_lambdas*(params.tmp_param.k- params.tmp_param.fixcluster0))+(params.tmp_param.k-1);

					if( cafe_param->parameters ) memory_free(cafe_param->parameters);
					cafe_param->parameters = NULL;
					cafe_param->parameters = (double*)memory_new(cafe_param->num_params, sizeof(double));
					if (cafe_param->k_weights) { memory_free(cafe_param->k_weights);}
					cafe_param->k_weights = NULL;
					cafe_param->k_weights = (double*) memory_new(cafe_param->k, sizeof(double));

					//cafe_param->lambda = cafe_param->parameters;
				}
				else {	// search whole dataset branch specific
					cafe_param->num_params = cafe_param->num_lambdas;
					
					if( cafe_param->parameters ) memory_free(cafe_param->parameters);
					cafe_param->parameters = NULL;
					cafe_param->parameters = (double*)memory_new(cafe_param->num_params, sizeof(double));
					
					//cafe_param->lambda = cafe_param->parameters;
				}
			}
			else {
				cafe_param->num_lambdas = params.tmp_param.num_lambdas = 1;
				if (params.tmp_param.k > 0) {
					cafe_param->k = params.tmp_param.k;
					cafe_param->fixcluster0 = params.tmp_param.fixcluster0;
					cafe_param->num_params = (params.tmp_param.num_lambdas*(params.tmp_param.k- params.tmp_param.fixcluster0))+(params.tmp_param.k-1);
					
					if( cafe_param->parameters ) memory_free(cafe_param->parameters);
					cafe_param->parameters = NULL;
					cafe_param->parameters = (double*)memory_new(cafe_param->num_params, sizeof(double));
					if (cafe_param->k_weights) { memory_free(cafe_param->k_weights);}
					cafe_param->k_weights = NULL;
					cafe_param->k_weights = (double*) memory_new(cafe_param->k, sizeof(double));
					
					//cafe_param->lambda = cafe_param->parameters;
				}
				else {	// search whole dataset whole tree
					cafe_param->num_params = cafe_param->num_lambdas;
					
					if( cafe_param->parameters ) memory_free(cafe_param->parameters);
					cafe_param->parameters = NULL;
					cafe_param->parameters = (double*)memory_new(cafe_param->num_params, sizeof(double));
					
					//cafe_param->lambda = cafe_param->parameters;
				}
			}
			// search
			if (params.tmp_param.checkconv) { cafe_param->checkconv = 1; }
			if (params.each )
			{
				cafe_each_best_lambda_by_fminsearch(cafe_param,cafe_param->num_lambdas);
			}
			else
			{
				cafe_best_lambda_by_fminsearch(cafe_param, cafe_param->num_lambdas, cafe_param->k);
			}
            // scale back branch lengths by sum_branch_length.
            /*for( j = 0 ; j < ptree->nlist->size; j++ )
            {
                pPhylogenyNode pnode = (pPhylogenyNode)ptree->nlist->array[j];
                if ( pnode->branchlength > 0 )
                {
                    pnode->branchlength = pnode->branchlength*cafe_param->sum_branch_length;
                }
            }*/
		}
		else {
			if (params.tmp_param.lambda_tree != NULL) {
				// cafe_param->num_lambdas determined by lambda tree.
				if (params.tmp_param.k > 0) {	// search clustered branch specific
					cafe_param->k = params.tmp_param.k;
					cafe_param->fixcluster0 = params.tmp_param.fixcluster0;
					cafe_param->num_params = (params.tmp_param.num_lambdas*(params.tmp_param.k- params.tmp_param.fixcluster0))+(params.tmp_param.k-1);
					
					// check if the numbers of lambdas and proportions put in matches the number of parameters
					if (cafe_param->num_params != params.tmp_param.num_params) {
						fprintf( stderr, "ERROR(lambda): Number of parameters not correct. \n");
						fprintf( stderr, "the number of -l lambdas and -p proportions are %d they need to be %d\n", params.tmp_param.num_params, cafe_param->num_params );
						pString pstr = phylogeny_string(params.tmp_param.lambda_tree,NULL);
						fprintf( stderr, "based on the tree %s and -k clusters %d.\n", pstr->buf, cafe_param->k );
						string_free(pstr);
						return -1;						
					}
					
					// copy user input into parameters
					if( cafe_param->parameters ) memory_free(cafe_param->parameters);
					cafe_param->parameters = NULL;
					cafe_param->parameters = (double*)memory_new(cafe_param->num_params, sizeof(double));
					memcpy(cafe_param->parameters,params.tmp_param.lambda, sizeof(double)*params.tmp_param.num_lambdas*(params.tmp_param.k-params.tmp_param.fixcluster0));
					memcpy(&cafe_param->parameters[cafe_param->num_lambdas*(cafe_param->k-params.tmp_param.fixcluster0)], params.tmp_param.k_weights, sizeof(double)*(params.tmp_param.k-1));
					// prepare space for k_weights
					if ( cafe_param->k_weights ) memory_free(cafe_param->k_weights);
					cafe_param->k_weights = NULL;
					cafe_param->k_weights = (double*)memory_new(cafe_param->k-1, sizeof(double) );										
					
				}
				else {	// search whole dataset branch specific
					cafe_param->num_params = cafe_param->num_lambdas;
					
					// check if the numbers of lambdas and proportions put in matches the number of parameters
					if (cafe_param->num_params != params.tmp_param.num_params) {
						fprintf( stderr, "ERROR(lambda): Number of parameters not correct. \n");
						fprintf( stderr, "the number of -l lambdas are %d they need to be %d\n", params.tmp_param.num_params, cafe_param->num_params );
						pString pstr = phylogeny_string(params.tmp_param.lambda_tree,NULL);
						fprintf( stderr, "based on the tree %s \n", pstr->buf );
						string_free(pstr);
						return -1;						
					}
					
					// copy user input into parameters
					if( cafe_param->parameters ) memory_free(cafe_param->parameters);
					cafe_param->parameters = NULL;
					cafe_param->parameters = (double*)memory_new(cafe_param->num_params, sizeof(double));
					memcpy(cafe_param->parameters,params.tmp_param.lambda, sizeof(double)*params.tmp_param.num_lambdas);
				}
			}
			else {
				cafe_param->num_lambdas = params.tmp_param.num_lambdas = 1;	
				if (params.tmp_param.k > 0) {				// search clustered whole tree
					cafe_param->k = params.tmp_param.k;
					cafe_param->fixcluster0 = params.tmp_param.fixcluster0;
					cafe_param->num_params = (cafe_param->num_lambdas*(cafe_param->k-cafe_param->fixcluster0))+(cafe_param->k-1);
					
					// check if the numbers of lambdas and proportions put in matches the number of parameters
					if (cafe_param->num_params != params.tmp_param.num_params) {
						fprintf( stderr, "ERROR(lambda): Number of parameters not correct. \n");
						fprintf( stderr, "the number of -l lambdas and -p proportions are %d they need to be %d\n", params.tmp_param.num_params, cafe_param->num_params );
						fprintf( stderr, "based on the -k clusters %d.\n", cafe_param->k );
						return -1;						
					}
					
					// copy user input into parameters
					if( cafe_param->parameters ) memory_free(cafe_param->parameters);
					cafe_param->parameters = NULL;
					cafe_param->parameters = (double*)memory_new(cafe_param->num_params, sizeof(double));
					memcpy(cafe_param->parameters,params.tmp_param.lambda, sizeof(double)*params.tmp_param.num_lambdas*(params.tmp_param.k-params.tmp_param.fixcluster0));
					memcpy(&cafe_param->parameters[cafe_param->num_lambdas*(cafe_param->k-params.tmp_param.fixcluster0)], params.tmp_param.k_weights, sizeof(double)*(params.tmp_param.k-1));
					// prepare space for k_weights
					if ( cafe_param->k_weights ) memory_free(cafe_param->k_weights);
					cafe_param->k_weights = NULL;
					cafe_param->k_weights = (double*)memory_new(cafe_param->k-1, sizeof(double) );										
					
				}
				else {	// search whole dataset whole tree
					cafe_param->num_params = cafe_param->num_lambdas;
					
					// check if the numbers of lambdas and proportions put in matches the number of parameters
					if (cafe_param->num_params != params.tmp_param.num_params) {
						fprintf( stderr, "ERROR(lambda): Number of parameters not correct. \n");
						fprintf( stderr, "the number of -l lambdas are %d they need to be %d\n", params.tmp_param.num_params, cafe_param->num_params );
						return -1;						
					}
					
					// copy user input into parameters
					if( cafe_param->parameters ) memory_free(cafe_param->parameters);
					cafe_param->parameters = NULL;
					cafe_param->parameters = (double*)memory_new(cafe_param->num_params, sizeof(double));
					memcpy(cafe_param->parameters,params.tmp_param.lambda, sizeof(double)*params.tmp_param.num_lambdas);
				}
			}
			cafe_param->param_set_func(cafe_param, cafe_param->parameters);
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
		for ( i = 0 ; i < cafe_param->pfamily->flist->size ; i++ )
		{
			pCafeFamilyItem pitem = (pCafeFamilyItem)cafe_param->pfamily->flist->array[i];
			cafe_family_set_size(cafe_param->pfamily, i, pcafe);
			cafe_param->param_set_func(cafe_param,pitem->lambda);
			pString pstr = cafe_tree_string_with_familysize_lambda(pcafe);
			for ( j = 0 ; j < cafe_param->num_lambdas; j++ )
			{
				double a = pitem->lambda[j] * cafe_param->max_branch_length;
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
		if ( cafe_param->pfamily )
		{
			cafe_set_birthdeath_cache_thread(cafe_param);
		}
	}
	
	cafe_log(cafe_param,"DONE: Lambda Search or setting, for command:\n");

	std::ostringstream cmd;
	std::copy(tokens.begin(), tokens.end(),
		std::ostream_iterator<std::string>(cmd, " "));
	cafe_log(cafe_param,"%s\n", cmd.str().c_str());
	
	if (params.search && (cafe_param->k > 0)) {
		// print the cluster memberships
		cafe_family_print_cluster_membership(cafe_param);
	}
	}
	catch (const std::exception& ex)
	{
		fprintf(stderr, "%s\n", ex.what());
		return -1;
	}

	return 0;
}



