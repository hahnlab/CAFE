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
#include <iostream>
#include <numeric>

#include "lambda.h"
#include "cafe_commands.h"
#include "Globals.h"
#include "log_buffer.h"
#include "gene_family.h"

extern "C" {
	extern pCafeParam cafe_param;
	void cafe_shell_set_lambda(pCafeParam param, double* parameters);
	int __cafe_cmd_lambda_tree(pCafeParam param, char *arg1, char *arg2);
	double __cafe_best_lambda_search(double* plambda, void* args);
	double __cafe_cluster_lambda_search(double* parameters, void* args);
}

using namespace std;

bool is_out(Argument arg)
{
	return !strcmp(arg.opt.c_str(), "-o");
}

void lambda_arg_base::load(vector<Argument> pargs)
{
	for (size_t i = 0; i < pargs.size(); i++)
	{
		pArgument parg = &pargs[i];

		// Search for whole family 
		if (!strcmp(parg->opt.c_str(), "-s"))
		{
			search = true;
		}
		else if (!strcmp(parg->opt.c_str(), "-checkconv"))
		{
			checkconv = true;
		}
		else if (!strcmp(parg->opt.c_str(), "-t"))
		{
            char *arg1 = NULL, *arg2 = NULL;
            if (parg->argc == 2)
            {
                arg1 = parg->argv[0];
                arg2 = parg->argv[1];
            }
            else
            {
            }
			bdone = __cafe_cmd_lambda_tree(cafe_param, parg->argv[0], parg->argc > 1 ? parg->argv[1] : NULL);
			if (bdone < 0) {
				throw std::exception();
			}
			pString pstr = phylogeny_string(cafe_param->lambda_tree, NULL);
			cafe_log(cafe_param, "Lambda Tree: %s\n", pstr->buf);
			string_free(pstr);
			lambda_tree = cafe_param->lambda_tree;
			lambdas.resize(cafe_param->num_lambdas);

		}
		else if (!strcmp(parg->opt.c_str(), "-l"))
		{
			get_doubles_array(lambdas, parg);
			num_params += lambdas.size();
			lambda_type = MULTIPLE_LAMBDAS;

		}
		else if (!strcmp(parg->opt.c_str(), "-p"))
		{
			get_doubles_array(k_weights, parg);
			num_params += k_weights.size();
		}
		else if (!strcmp(parg->opt.c_str(), "-k"))
		{
			int k = 0;
			sscanf(parg->argv[0], "%d", &k);
			k_weights.resize(k);
		}
		else if (!strcmp(parg->opt.c_str(), "-f"))
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
		if (!strcmp(parg->opt.c_str(), "-v"))
		{
			sscanf(parg->argv[0], "%lf", &vlambda);
			lambda_type = SINGLE_LAMBDA;
		} 
		else if (!strcmp(parg->opt.c_str(), "-r"))
		{
			range.resize(parg->argc);
			int j;
			for (j = 0; j < parg->argc; j++)
			{
				sscanf(parg->argv[j], "%lf:%lf:%lf", &range[j].start, &range[j].step, &range[j].end);
			}
		}
		else if (!strcmp(parg->opt.c_str(), "-e"))
		{
			write_files = true;
			each = true;
		}
		else if (!strcmp(parg->opt.c_str(), "-o"))
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
		input_values_construct(&param->input, param->num_params);
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
	param->num_params = params.get_num_params();

	params.validate_parameter_count(param->num_params);

	initialize_params_and_k_weights(param, INIT_PARAMS | INIT_KWEIGHTS);

	input_values_set_lambdas(&param->input, &params.lambdas[0], params.lambdas.size());

	int num_k_weights = params.k_weights.size() - 1;
	int first_k_weight = (param->num_lambdas*(param->parameterized_k_value - params.fixcluster0));
	input_values_set_k_weights(&param->input, &params.k_weights[0], first_k_weight, num_k_weights);

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
			input_values_set_lambdas(&param->input, &params.lambdas[0], params.lambdas.size());
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
			input_values_set_lambdas(&param->input, &params.lambdas[0], param->num_lambdas);
		}
	}

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

  log_buffer buf(&globals.param);
  ostream ost(&buf);
  
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

        std::vector<double> prior_rfsize;
		cafe_set_prior_rfsize_empirical(param, prior_rfsize);
        globals.param.prior_rfsize = (double *)memory_new(prior_rfsize.size(), sizeof(double));
        std::copy(prior_rfsize.begin(), prior_rfsize.end(), globals.param.prior_rfsize);

		param->num_params = param->num_lambdas;
        
		initialize_params_and_k_weights(param, INIT_PARAMS);

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
    std::vector<double> prior_rfsize;
    cafe_set_prior_rfsize_empirical(param, prior_rfsize);
    globals.param.prior_rfsize = (double *)memory_new(prior_rfsize.size(), sizeof(double));
    std::copy(prior_rfsize.begin(), prior_rfsize.end(), globals.param.prior_rfsize);

	// search or set
	if (params.search) {
		lambda_search(param, params);
	}
	else {
		lambda_set(param, params);
		cafe_shell_set_lambda(param, param->input.parameters);
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
			cafe_family_set_size(param->pfamily, pitem, pcafe);
            cafe_shell_set_lambdas(param,pitem->lambda);
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
		log_cluster_membership(param->pfamily, param->parameterized_k_value, param->p_z_membership, ost);
	}

	return 0;
}

double* cafe_best_lambda_by_fminsearch(pCafeParam param, int lambda_len, int k)
{
	int i, j;
	int max_runs = 10;
	double* scores = (double *)memory_new(max_runs, sizeof(double));
	int converged = 0;
	int runs = 0;

	do
	{

		if (param->num_params > 0)
		{
			int kfix = k - param->fixcluster0;
			double max_branch_length = param->max_branch_length;
			double *k_weights = param->k_weights;
			input_values_randomize(&param->input, param->num_lambdas, param->num_mus, param->parameterized_k_value,
				kfix, max_branch_length, k_weights);
		}

		copy_range_to_tree(param->pcafe, &param->family_size);

		pFMinSearch pfm;
		if (k > 0) {
			pfm = fminsearch_new_with_eq(__cafe_cluster_lambda_search, param->num_params, param);
			pfm->tolx = 1e-5;
			pfm->tolf = 1e-5;
		}
		else {
			pfm = fminsearch_new_with_eq(__cafe_best_lambda_search, lambda_len, param);
			pfm->tolx = 1e-6;
			pfm->tolf = 1e-6;
		}

        // pass a copy of input_parameters to fminsearch to guarantee that parameters won't be changed inside 
        // the search function
        std::vector<double> starting_values(param->num_params);
        std::copy(param->input.parameters, param->input.parameters + param->num_params, starting_values.begin());

        fminsearch_min(pfm, &starting_values[0]);

        double *re = fminsearch_get_minX(pfm);
		for (i = 0; i < param->num_params; i++) param->input.parameters[i] = re[i];

		double current_p = param->input.parameters[(lambda_len)*(k - param->fixcluster0)];
		double prev_p;
		if (k>0)
		{
			do {
				double* sumofweights = (double*)memory_new(param->parameterized_k_value, sizeof(double));
				for (i = 0; i < param->pfamily->flist->size; i++) {
					for (j = 0; j<k; j++) {
						sumofweights[j] += param->p_z_membership[i][j];
					}
				}
				for (j = 0; j<k - 1; j++) {
					param->input.parameters[(lambda_len)*(k - param->fixcluster0) + j] = sumofweights[j] / param->pfamily->flist->size;
				}
				memory_free(sumofweights);

				fminsearch_min(pfm, param->input.parameters);

				double *re = fminsearch_get_minX(pfm);
				for (i = 0; i < param->num_params; i++) param->input.parameters[i] = re[i];

				prev_p = current_p;
				current_p = param->input.parameters[(lambda_len)*(k - param->fixcluster0)];
			} while (current_p - prev_p > pfm->tolx);
		}

		cafe_log(param, "\n");
		cafe_log(param, "Lambda Search Result: %d\n", pfm->iters);
		if (k > 0) {
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			if (param->fixcluster0) {
				strncat(buf, "0,", 2);
				string_pchar_join_double(buf, ",", param->num_lambdas*(param->parameterized_k_value - param->fixcluster0), param->input.parameters);
			}
			else {
				string_pchar_join_double(buf, ",", param->num_lambdas*param->parameterized_k_value, param->input.parameters);
			}
			cafe_log(param, "Lambda : %s\n", buf);
			buf[0] = '\0';
			if (param->parameterized_k_value > 0) {
				string_pchar_join_double(buf, ",", param->parameterized_k_value, param->k_weights);
				cafe_log(param, "p : %s\n", buf);
				cafe_log(param, "p0 : %f\n", param->input.parameters[param->num_lambdas*(param->parameterized_k_value - param->fixcluster0) + 0]);
			}
			cafe_log(param, "Score: %f\n", *pfm->fv);
		}
		else {
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			string_pchar_join_double(buf, ",", param->num_lambdas, param->input.parameters);
			cafe_log(param, "Lambda : %s & Score: %f\n", buf, *pfm->fv);
		}
		if (runs > 0) {
			double minscore = __min(scores, runs);
			if (abs(minscore - (*pfm->fv)) < 10 * pfm->tolf) {
				converged = 1;
			}
		}
		scores[runs] = *pfm->fv;
		fminsearch_free(pfm);

		copy_range_to_tree(param->pcafe, &param->family_size);

		runs++;

	} while (param->checkconv && !converged && runs<max_runs);

	//	string_free(pstr);
	if (param->checkconv) {
		if (converged) {
			cafe_log(param, "score converged in %d runs.\n", runs);
		}
		else {
			cafe_log(param, "score failed to converge in %d runs.\n", max_runs);
		}
	}
	memory_free(scores);
	return param->input.parameters;
}

posterior compute_posterior(pCafeFamilyItem pitem, pCafeTree pcafe, const std::vector<double>& prior_rfsize)
{
    posterior result;
    compute_tree_likelihoods(pcafe);

    double *likelihood = get_likelihoods(pcafe);		// likelihood of the whole tree = multiplication of likelihood of all nodes

    result.max_likelihood = __max(likelihood, pcafe->rfsize);			// this part find root size condition with maxlikelihood for each family			
    if (pitem->maxlh < 0)
    {
        pitem->maxlh = __maxidx(likelihood, pcafe->rfsize);
    }
    // get posterior by adding lnPrior to lnLikelihood
    std::vector<double> posterior(pcafe->rfsize);
    for (int j = 0; j < pcafe->rfsize; j++)	// j: root family size
    {
        // likelihood and posterior both starts from 1 instead of 0 
        posterior[j] = exp(log(likelihood[j]) + log(prior_rfsize[j]));	//prior_rfsize also starts from 1
    }

    // this part find root size condition with maxlikelihood for each family			
    result.max_posterior = *std::max_element(posterior.begin(), posterior.end());

    return result;
}

double get_posterior(pCafeFamily pfamily, pCafeTree pcafe, family_size_range*range, std::vector<double>& ML, std::vector<double>& MAP, std::vector<double>& prior_rfsize, int quiet)
{
	int i;
	double score = 0;
	for (i = 0; i < pfamily->flist->size; i++)	// i: family index
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)pfamily->flist->array[i];

        if (pitem->ref < 0 || pitem->ref == i)
		{
			cafe_family_set_size(pfamily, pitem, pcafe);	// this part is just setting the leave counts.
			posterior p = compute_posterior(pitem, pcafe, prior_rfsize);
            ML[i] = p.max_likelihood;
            MAP[i] = p.max_posterior;

		}
		else
		{
            ML[i] = ML[pitem->ref];
			MAP[i] = MAP[pitem->ref];
		}
		if (ML[i] == 0)
		{
			ostringstream ost;
			ost << "WARNING: Calculated posterior probability for family " << pitem->id << " = 0" << endl;
			throw std::runtime_error(ost.str());
		}
		score += log(MAP[i]);			// add log-posterior across all families
	}
	return score;
}

double __cafe_best_lambda_search(double* plambda, void* args)
{
	int i;
	pCafeParam param = (pCafeParam)args;
	pCafeTree pcafe = (pCafeTree)param->pcafe;
	double score = 0;
	int skip = 0;
	for (i = 0; i < param->num_lambdas; i++)
	{
		if (plambda[i] < 0)
		{
			skip = 1;
			score = log(0);
			break;
		}
	}
	if (!skip)
	{
        cafe_shell_set_lambdas(param, plambda);

		reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, &param->family_size);
		try
		{
            int num_families = param->pfamily->flist->size;
            std::vector<double> ml(num_families), map(num_families), pr(FAMILYSIZEMAX);
            copy(param->ML, param->ML + num_families, ml.begin());
            copy(param->MAP, param->MAP + num_families, map.begin());
            copy(param->prior_rfsize, param->prior_rfsize + FAMILYSIZEMAX, pr.begin());
            score = get_posterior(param->pfamily, param->pcafe, &param->family_size, ml, map, pr, param->quiet);
		}
		catch (std::runtime_error& e)
		{
            if (!param->quiet)
            {
                std::cerr << e.what();
            }
			score = log(0);
		}
		cafe_free_birthdeath_cache(pcafe);
	}
	char buf[STRING_STEP_SIZE];
	buf[0] = '\0';
	string_pchar_join_double(buf, ",", param->num_lambdas, plambda);
	cafe_log(param, "Lambda : %s & Score: %f\n", buf, score);
	cafe_log(param, ".");
	return -score;
}

double __lnLPoisson(double* plambda, void* data)
{
    double score = 0;
    double lambda = plambda[0];
    std::vector<int>* leaf_sizes = (std::vector<int> *)data;

    for (size_t i = 0; i<leaf_sizes->size(); i++) {
        int x = leaf_sizes->at(i);
        double ll = poisspdf((double)x, lambda);
        if (std::isnan(ll)) {
            ll = 0;
        }
        score += log(ll);
    }

    return -score;
}

std::vector<int> collect_leaf_sizes(pCafeFamily pfamily)
{
    // estimate the distribution of family size based on observed leaf counts.
    // first collect all leaves sizes into an ArrayList.
    std::vector<int> pLeavesSize;
    for (int idx = 0; idx<pfamily->flist->size; idx++) {
        pCafeFamilyItem pitem = (pCafeFamilyItem)pfamily->flist->array[idx];
        for (int i = 0; i < pfamily->num_species; i++)
        {
            if (pfamily->index[i] < 0) continue;
            if (pitem->count[i] > 0) {		// ignore the zero counts ( we condition that rootsize is at least one )
                pLeavesSize.push_back(pitem->count[i]-1);
            }
        }
    }

    return pLeavesSize;
}

poisson_lambda find_poisson_lambda(pCafeFamily pfamily)
{
    poisson_lambda result;
    int i = 0;

    auto leaf_sizes = collect_leaf_sizes(pfamily);

    // now estimate parameter based on data and distribution (poisson or gamma). 
    pFMinSearch pfm;
    int num_params = 1;
    pfm = fminsearch_new_with_eq(__lnLPoisson, num_params, &leaf_sizes);
    //int num_params = 2;
    //pfm = fminsearch_new_with_eq(__lnLGamma,num_params,pLeavesSize);
    pfm->tolx = 1e-6;
    pfm->tolf = 1e-6;
    double* parameters = (double *)memory_new(num_params, sizeof(double));
    for (i = 0; i < num_params; i++) parameters[i] = unifrnd();
    fminsearch_min(pfm, parameters);
    double *re = fminsearch_get_minX(pfm);
    for (i = 0; i < num_params; i++) parameters[i] = re[i];

    result.parameters = parameters;
    result.num_params = num_params;
    result.num_iterations = pfm->iters;
    result.score = *pfm->fv;

    // clean
    fminsearch_free(pfm);

    return result;
}


void cafe_set_prior_rfsize_poisson_lambda(std::vector<double>& prior_rfsize, int shift, double* lambda)
{
    // shift = param->pcafe->rootfamilysizes[0]
    int i;
    // calculate the prior probability for a range of root sizes.
    prior_rfsize.resize(FAMILYSIZEMAX);
    for (i = 0; i<FAMILYSIZEMAX; i++) {
        //param->prior_rfsize[i] = poisspdf(param->pcafe->rootfamilysizes[0]+i, parameters[0]);					// poisson
        prior_rfsize[i] = poisspdf(shift - 1 + i, lambda[0]);					// shifted poisson
        //param->prior_rfsize[i] = gampdf(param->pcafe->rootfamilysizes[0]+i, parameters[0], parameters[1]);	// gamma
    }
}


/// set empirical prior on rootsize based on the assumption that rootsize follows leaf size distribution
double cafe_set_prior_rfsize_empirical(pCafeParam param, std::vector<double>& prior_rfsize)
{
    poisson_lambda result = find_poisson_lambda(param->pfamily);
    cafe_log(param, "Empirical Prior Estimation Result: (%d iterations)\n", result.num_iterations);
    cafe_log(param, "Poisson lambda: %f & Score: %f\n", result.parameters[0], result.score);

    double *prior_poisson_lambda = (double *)memory_new_with_init(result.num_params, sizeof(double), (void*)result.parameters);
    //cafe_log(param,"Gamma alpha: %f, beta: %f & Score: %f\n", parameters[0], parameters[1], *pfm->fv);	

    // set rfsize based on estimated prior
    cafe_set_prior_rfsize_poisson_lambda(prior_rfsize, param->pcafe->rootfamilysizes[0], prior_poisson_lambda);

    memory_free(result.parameters);
    return 0;
}

