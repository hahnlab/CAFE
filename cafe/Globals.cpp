#include <stdexcept>

#include "Globals.h"
#include "viterbi.h"

extern "C" {
#include "family.h"
#include "cafe.h"

	void cafe_shell_set_lambda(pCafeParam param, double* parameters);
	extern pCafeParam cafe_param;
	extern pBirthDeathCacheArray probability_cache;

}
/**
* \brief Initializes the global \ref cafe_param that holds the data acted upon by cafe. Called at program startup.
*
*/
Globals::Globals() : viterbi(new viterbi_parameters())
{
	param.family_size.root_min = 1;
	param.family_size.root_max = 1;
	param.family_size.min = 0;
	param.family_size.max = 1;
	param.param_set_func = cafe_shell_set_lambda;
	param.flog = stdout;
	param.num_threads = 1;
	param.num_random_samples = 1000;
	param.pvalue = 0.01;
	param.quiet = 0;
	param.pcafe = NULL;
	param.pfamily = NULL;
	param.lambda = NULL;
	param.lambda_tree = NULL;
	param.ML = NULL;
	param.MAP = NULL;
	param.prior_rfsize = NULL;
	param.parameters = NULL;
	param.str_fdata = NULL;
	viterbi->viterbiPvalues = NULL;
	viterbi->cutPvalues = NULL;
	mu_tree = NULL;

	cafe_param = &param;
}

Globals::~Globals()
{
	delete viterbi;
}

void Globals::Clear(int btree_skip)
{
	if (param.str_fdata)
	{
		string_free(param.str_fdata);
		param.str_fdata = NULL;
	}
	if (param.ML)
	{
		memory_free(param.ML);
		param.ML = NULL;
	}
	if (param.MAP)
	{
		memory_free(param.MAP);
		param.MAP = NULL;
	}
	if (param.prior_rfsize)
	{
		memory_free(param.prior_rfsize);
		param.prior_rfsize = NULL;
	}

	int nnodes = param.pcafe ? ((pTree)param.pcafe)->nlist->size : 0;
	viterbi_parameters_clear(viterbi, nnodes);
	if (!btree_skip && param.pcafe)
	{
		if (probability_cache)
		{
			cafe_free_birthdeath_cache(param.pcafe);
		}
		cafe_tree_free(param.pcafe);
		//memory_free( param.branchlengths_sorted );
		//param.branchlengths_sorted = NULL;
		memory_free(param.old_branchlength);
		param.old_branchlength = NULL;
		param.pcafe = NULL;
	}
	if (param.pfamily)
	{
		cafe_family_free(param.pfamily);
		param.pfamily = NULL;
	}
	if (param.parameters)
	{
		memory_free(param.parameters);
		param.parameters = NULL;
	}
	if (param.lambda)
	{
		//memory_free( param.lambda ); param.lambda points to param.parameters
		param.lambda = NULL;
	}
	if (param.mu)
	{
		//memory_free( param.mu );	param.mu points to param.parameters
		param.mu = NULL;
	}
	if (param.lambda_tree)
	{
		phylogeny_free(param.lambda_tree);
		param.lambda_tree = NULL;
	}
	if (mu_tree)
	{
		phylogeny_free(mu_tree);
		mu_tree = NULL;
	}

	param.eqbg = 0;
	param.posterior = 0;
	param.num_params = 0;
	param.num_lambdas = 0;
	param.num_mus = 0;
	param.parameterized_k_value = 0;
	param.fixcluster0 = 0;
	param.family_size.root_min = 0;
	param.family_size.root_max = 1;
	param.family_size.min = 0;
	param.family_size.max = 1;
	param.param_set_func = cafe_shell_set_lambda;
	param.num_threads = 1;
	param.num_random_samples = 1000;
	param.pvalue = 0.01;
}

void Globals::Prepare()
{
	if (param.pfamily == NULL || param.pcafe == NULL)
	{
		throw std::runtime_error("ERROR(lambda): Please load family (\"load\") and cafe tree (\"tree\") before running \"lambda\" command.");
	}

	param.lambda = NULL;
	param.mu = NULL;
	if (param.lambda_tree) {
		phylogeny_free(param.lambda_tree);
		param.lambda_tree = NULL;
	}
	if (mu_tree)
	{
		phylogeny_free(mu_tree);
		mu_tree = NULL;
	}
	param.num_lambdas = -1;
	param.num_mus = -1;
	param.parameterized_k_value = 0;
	param.param_set_func = cafe_shell_set_lambda;

}


