#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <cmath>
#include <iostream>

#include "lambdamu.h"
#include "Globals.h"
#include "cafe_commands.h"
#include "lambda.h"
#include "log_buffer.h"
#include "gene_family.h"

using namespace std;

extern "C"
{
	extern pCafeParam cafe_param;
	void cafe_shell_set_lambda_mu(pCafeParam param, double* parameters);

#include "cafe.h"
}

void lambdamu_args::load(std::vector<Argument> args)
{
	lambda_arg_base::load(args);

	for (size_t i = 0; i < args.size(); i++)
	{
		pArgument parg = &args[i];

		if (!strcmp(parg->opt.c_str(), "-t"))
		{
			cafe_param->num_mus = cafe_param->num_lambdas;
			mus.resize(cafe_param->num_mus);
		}
		else if (!strcmp(parg->opt.c_str(), "-m"))
		{
			get_doubles_array(mus, parg);
			num_params += mus.size();
		}
		else if (!strcmp(parg->opt.c_str(), "-eqbg"))
		{
			eqbg = 1;
		}
	}
}

int lambdamu_args::get_num_params() const
{
	if (lambda_tree == NULL)
	{
		if (k_weights.size() > 0) {
			return (lambdas.size()*(k_weights.size() - fixcluster0)) +
				(mus.size()*(k_weights.size() - fixcluster0)) +
				(k_weights.size() - 1);
		}
		else
			return lambdas.size() + (mus.size() - eqbg);
	}
	else
	{
		if (k_weights.size() > 0) {

			return (lambdas.size()*(k_weights.size() - fixcluster0)) +
				((mus.size() - eqbg)*(k_weights.size() - fixcluster0)) +
				(k_weights.size() - 1);
		}
		else
		{
			return lambdas.size() + mus.size();
		}
	}
}


void lambdamu_prepare_search(pCafeParam param, lambdamu_args& params)
{
	// prepare parameters
	if (params.lambda_tree != NULL) {
		// param->num_lambdas determined by lambda tree.
		param->eqbg = params.eqbg;
	}
	else {
		param->num_lambdas = 1;
		params.lambdas.resize(1);
		param->num_mus = 1;
		params.mus.resize(1);
		if (params.eqbg) {
			throw runtime_error("ERROR(lambdamu): Cannot use option eqbg without specifying a lambda tree. \n");
		}
	}

	if (params.k_weights.size() > 0) {
		param->parameterized_k_value = params.k_weights.size();
		param->fixcluster0 = params.fixcluster0;
	}

	param->num_params = params.get_num_params();
	initialize_params_and_k_weights(param, INIT_PARAMS | (params.k_weights.empty() ? 0 : INIT_KWEIGHTS));
	if (params.checkconv) { param->checkconv = 1; }
}

void lambda_arg_base::validate_parameter_count(int expected)
{
	ostringstream ost;
	// check if the numbers of lambdas and proportions put in matches the number of parameters
	if (expected != num_params)
	{
		ost << "ERROR (" << command() << "): The total number of parameters was not correct.\n";
		ost << "The total number of " << args();
		if (k_weights.size() > 0)
			ost << " and proportions (-p)";
		ost << " are " << num_params << " but " << expected << " were expected\n";
		ost << "based on the ";
		if (lambda_tree)
		{
			pString pstr = phylogeny_string(lambda_tree, NULL);
			ost << "tree " << pstr->buf;
			string_free(pstr);
		}
		if (k_weights.size() > 0)
		{
			if (lambda_tree)
				ost << " and ";
			ost << "the " << k_weights.size() << " clusters (-k).";
		}
		ost << "\n";
		throw runtime_error(ost.str());
	}
}

void lambdamu_set(pCafeParam param, lambdamu_args& params)
{
	if (params.lambda_tree != NULL) {
		// param->num_lambdas determined by lambda tree.
		param->eqbg = params.eqbg;
		if (!params.k_weights.empty()) {	// search clustered branch specific
			param->parameterized_k_value = params.k_weights.size();
			param->fixcluster0 = params.fixcluster0;
			param->num_params = params.get_num_params();

			params.validate_parameter_count(param->num_params);

			// copy user input into parameters
			initialize_params_and_k_weights(param, INIT_PARAMS | INIT_KWEIGHTS);

			// Not sure what values of fixcluster0 might cause this to NOT overrun
			size_t num_lambdas = params.lambdas.size()*(params.k_weights.size() - params.fixcluster0);
			size_t num_mus = ((params.mus.size() - params.eqbg)*(params.k_weights.size() - params.fixcluster0));

			input_values_set_lambdas(&param->input, &params.lambdas[0], min(num_lambdas, params.lambdas.size()));

			int first_mu = param->num_lambdas*(params.k_weights.size() - params.fixcluster0);
			input_values_set_mus(&param->input, &params.mus[0], first_mu, min(num_mus, params.mus.size()));

			int first_k_weight = (param->num_lambdas*(params.k_weights.size() - params.fixcluster0)) +
				((params.mus.size() - params.eqbg)*(params.k_weights.size() - params.fixcluster0));

			input_values_set_mus(&param->input, &params.k_weights[0], first_k_weight, min(num_mus, params.k_weights.size()-1));
		}
		else {	// search whole dataset branch specific
			param->num_params = params.lambdas.size() + (params.mus.size() - params.eqbg);

			params.validate_parameter_count(param->num_params);

			// copy user input into parameters
			initialize_params_and_k_weights(param, INIT_PARAMS);

			input_values_set_lambdas(&param->input, &params.lambdas[0], params.lambdas.size());
			input_values_set_mus(&param->input, &params.mus[0], param->num_lambdas, params.mus.size() - params.eqbg);
		}
	}
	else {
		param->num_lambdas = 1;
		params.lambdas.resize(1);
		param->num_mus = 1;
		params.mus.resize(1);
		if (params.eqbg) {
			throw runtime_error("ERROR(lambdamu): Cannot use option eqbg without specifying a lambda tree. \n");
		}
		if (params.k_weights.size() > 0) {				// search clustered whole tree
			param->parameterized_k_value = params.k_weights.size();
			param->fixcluster0 = params.fixcluster0;
			param->num_params = params.get_num_params();

			params.validate_parameter_count(param->num_params);

			// copy user input into parameters
			initialize_params_and_k_weights(param, INIT_PARAMS | INIT_KWEIGHTS);

			size_t num_lambdas = params.lambdas.size()*(params.k_weights.size() - params.fixcluster0);
			copy(params.lambdas.begin(), params.lambdas.begin() + min(params.lambdas.size(), num_lambdas), param->input.parameters);

			int first_mu = param->num_lambdas*(params.k_weights.size() - params.fixcluster0);
			size_t num_mus = params.mus.size()*(params.k_weights.size() - params.fixcluster0);
			copy(params.mus.begin(), params.mus.begin() + min(params.mus.size(), num_mus), param->input.parameters + first_mu);

			int first_k_weight = param->num_lambdas*(params.k_weights.size() - params.fixcluster0) + params.mus.size()*(param->parameterized_k_value - params.fixcluster0);
			copy(params.k_weights.begin(), params.k_weights.end() - 1, param->input.parameters + first_k_weight);
		}
		else {	// search whole dataset whole tree
			param->num_params = params.lambdas.size() + params.mus.size();

			// check if the numbers of lambdas and proportions put in matches the number of parameters
			params.validate_parameter_count(param->num_params);

			// copy user input into parameters
			initialize_params_and_k_weights(param, INIT_PARAMS);

			copy(params.lambdas.begin(), params.lambdas.end(), param->input.parameters);
			copy(params.mus.begin(), params.mus.end(), param->input.parameters + params.lambdas.size());
		}
	}
}

int cafe_cmd_lambdamu(Globals& globals, std::vector<std::string> tokens)
{
	pCafeParam param = &globals.param;
	log_buffer buf(&globals.param);
	ostream log(&buf);

	prereqs(param, REQUIRES_FAMILY | REQUIRES_TREE);

	vector<Argument> args = build_argument_list(tokens);
	globals.Prepare();

	param->optimizer_init_type = LAMBDA_MU;

	lambdamu_args params;
	params.load(args);

	// copy parameters collected to param based on the combination of options.
	param->posterior = 1;
	// set rootsize prior based on leaf size
    std::vector<double> prior_rfsize;
    cafe_set_prior_rfsize_empirical(param, prior_rfsize);
    globals.param.prior_rfsize = (double *)memory_new(prior_rfsize.size(), sizeof(double));
    std::copy(prior_rfsize.begin(), prior_rfsize.end(), globals.param.prior_rfsize);

	// search or set
	if (params.search) {
		lambdamu_prepare_search(param, params);
		// search
		best_lambda_mu_by_fminsearch(param, param->num_lambdas, param->num_mus, param->parameterized_k_value, log);

	}
	else {
		lambdamu_set(param, params);
		cafe_shell_set_lambda_mu(param, param->input.parameters);
	}

	if (param->pfamily)
	{
		reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, &param->family_size);
	}

  log_complete(param, tokens, params.search, true);

  return 0;
}

double cafe_cluster_lambda_mu_search(double* parameters, void* args)
{
	int i;
	pCafeParam param = (pCafeParam)args;
	pCafeTree pcafe = (pCafeTree)param->pcafe;
	double score = 0;
	int skip = 0;
	for (i = 0; i < param->num_params; i++)
	{
		if (parameters[i] < 0)
		{
			skip = 1;
			score = log(0);
			break;
		}
	}
	if (!skip)
	{
        cafe_shell_set_lambdas(param, parameters);

		reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, &param->family_size);
		score = cafe_get_clustered_posterior(param, param->ML, param->MAP, param->prior_rfsize);
		cafe_free_birthdeath_cache(pcafe);
		cafe_tree_node_free_clustered_likelihoods(param);
	}
	char buf[STRING_STEP_SIZE];
	buf[0] = '\0';
	for (i = 0; i<param->num_lambdas; i++) {
		string_pchar_join_double(buf, ",", (param->parameterized_k_value - param->fixcluster0), &parameters[i*(param->parameterized_k_value - param->fixcluster0)]);
		fprintf(stdout, "Lambda branch %d: %s\n", i, buf);
		buf[0] = '\0';
	}
	for (i = 0; i<param->num_mus; i++) {
		string_pchar_join_double(buf, ",", (param->parameterized_k_value - param->fixcluster0), &parameters[param->num_lambdas*(param->parameterized_k_value - param->fixcluster0) + i*(param->parameterized_k_value - param->fixcluster0)]);
		fprintf(stdout, "Mu branch %d: %s \n", i, buf);
		buf[0] = '\0';
	}
	if (param->parameterized_k_value > 0) {
		string_pchar_join_double(buf, ",", param->parameterized_k_value, param->k_weights);
		fprintf(stdout, "p : %s\n", buf);
	}
	//cafe_log(param, "Score: %f\n", score);
	fprintf(stdout, ".");
	return -score;
}



// need to define a function with lambda and mu. 
// this function is provided as the equation to fmin search.
// also need to make a new param_set_func that includes mu.

double cafe_best_lambda_mu_search(double* parameters, void* args)
{
	int i;
	pCafeParam param = (pCafeParam)args;
	pCafeTree pcafe = (pCafeTree)param->pcafe;
	double score = 0;
	int skip = 0;
	for (i = 0; i < param->num_params; i++	)
	{
		if (parameters[i] < 0)
		{
			skip = 1;
			score = log(0);
			break;
		}
	}
	if (!skip)
	{
        cafe_shell_set_lambdas(param, parameters);

		reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, &param->family_size);
		try
		{
            std::vector<double> pr(FAMILYSIZEMAX);
            copy(param->prior_rfsize, param->prior_rfsize + FAMILYSIZEMAX, pr.begin());
            score = get_posterior(param->pfamily, param->pcafe, pr);
		}
		catch (std::runtime_error& e)
		{
			std::cerr << e.what() << endl;
			score = log(0);
		}

		cafe_free_birthdeath_cache(pcafe);
	}
	char buf[STRING_STEP_SIZE];
	buf[0] = '\0';
	string_pchar_join_double(buf, ",", param->num_lambdas, parameters);
	cafe_log(param, "Lambda : %s ", buf, score);
	buf[0] = '\0';
	string_pchar_join_double(buf, ",", param->num_mus - param->eqbg, parameters + param->num_lambdas);
	cafe_log(param, "Mu : %s & Score: %f\n", buf, score);
	cafe_log(param, ".");
	return -score;
}


void best_lambda_mu_by_fminsearch(pCafeParam param, int lambda_len, int mu_len, int k, ostream& log)
{
	int i;
	int max_runs = 10;
	vector<double> scores(max_runs);
	int converged = 0;
	int runs = 0;

	do
	{
		if (param->num_params > 0)
		{
			int kfix = k - param->fixcluster0;
			double *k_weights = param->k_weights;

			input_values_randomize(&param->input, param->num_lambdas, param->num_mus, param->parameterized_k_value,
				kfix, max_branch_length((pTree)param->pcafe), k_weights);
		}

		copy_range_to_tree(param->pcafe, &param->family_size);

		pFMinSearch pfm;
		if (k > 0) {
			pfm = fminsearch_new_with_eq(cafe_cluster_lambda_mu_search, param->num_params, param);
		}
		else {
			pfm = fminsearch_new_with_eq(cafe_best_lambda_mu_search, param->num_params, param);
		}
		pfm->tolx = 1e-6;
		pfm->tolf = 1e-6;
		fminsearch_min(pfm, param->input.parameters);
		double *re = fminsearch_get_minX(pfm);
		for (i = 0; i < param->num_params; i++) param->input.parameters[i] = re[i];

		log << "\n";
		log << "Lambda Search Result: " << pfm->iters << "\n";
		// print
		if (k>0) {
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			for (i = 0; i<param->num_lambdas; i++) {
				if (param->fixcluster0) {
					strncat(buf, "0,", 2);
					string_pchar_join_double(buf, ",", (param->parameterized_k_value - param->fixcluster0), &param->input.parameters[i*(param->parameterized_k_value - param->fixcluster0)]);
				}
				else {
					string_pchar_join_double(buf, ",", param->parameterized_k_value, &param->input.parameters[i*param->parameterized_k_value]);
				}
				log << "Lambda branch " << i << ": " << buf << "\n";
				buf[0] = '\0';
			}
			for (i = 0; i<param->num_mus - param->eqbg; i++) {
				if (param->fixcluster0) {
					strncat(buf, "0,", 2);
					string_pchar_join_double(buf, ",", (param->parameterized_k_value - param->fixcluster0), &param->input.parameters[param->num_lambdas*(param->parameterized_k_value - param->fixcluster0) + i*(param->parameterized_k_value - param->fixcluster0)]);
				}
				else {
					string_pchar_join_double(buf, ",", param->parameterized_k_value, &param->input.parameters[param->num_lambdas*param->parameterized_k_value + i*param->parameterized_k_value]);
				}
				log << "Mu branch " << i << ": " << buf << "\n";
				buf[0] = '\0';
			}
			if (param->parameterized_k_value > 0) {
				string_pchar_join_double(buf, ",", param->parameterized_k_value, param->k_weights);
				log << "p : " << buf << "\n";
				log << "p0 : " << param->input.parameters[param->num_lambdas*(param->parameterized_k_value - param->fixcluster0) + (param->num_mus - param->eqbg)*(param->parameterized_k_value - param->fixcluster0) + 0] << "\n";
			}
			log << "Score: " << *pfm->fv << "\n";
		}
		else {
			char buf[STRING_STEP_SIZE];
			buf[0] = '\0';
			string_pchar_join_double(buf, ",", param->num_lambdas, param->input.parameters);
			log << "Lambda : " << buf << " & Score: " << *pfm->fv;
			buf[0] = '\0';
			string_pchar_join_double(buf, ",", param->num_mus - param->eqbg, param->input.parameters + param->num_lambdas);
			log << "Mu : " << buf << " & Score: " << *pfm->fv << "\n";
		}
		if (runs > 0) {
			double minscore = *std::min_element(scores.begin(), scores.end());
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
			log << "score converged in " << runs << " runs.\n";
		}
		else {
			log << "score failed to converge in " << max_runs << " runs.\n";
		}
	}
}

