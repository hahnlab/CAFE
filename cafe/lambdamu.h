#ifndef LAMBDAMU_H_3DE9A984_42A1_40DA_A5EE_877793436272
#define LAMBDAMU_H_3DE9A984_42A1_40DA_A5EE_877793436272

#include "lambda.h"

struct lambdamu_args : lambda_arg_base
{
	lambdamu_args() : lambda_arg_base(), eqbg(false) {}

	virtual void load(std::vector<Argument> pargs);

	std::vector<double> mus;
	bool eqbg;

	virtual const char* command() { return "lambdamu"; }
	virtual const char* args() { return "lambdas (-l) and mus (-m)"; }

	virtual int get_num_params() const;
};

void lambdamu_set(pCafeParam param, lambdamu_args& params);
void best_lambda_mu_by_fminsearch(pCafeParam param, int lambda_len, int mu_len, int k, std::ostream& log);
double cafe_cluster_lambda_mu_search(double* parameters, void* args);
double cafe_best_lambda_mu_search(double* parameters, void* args);

#endif
