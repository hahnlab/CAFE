#ifndef LAMBDA_H_5CC2A46B_4881_4C8F_BDE7_BE29023F2D15
#define LAMBDA_H_5CC2A46B_4881_4C8F_BDE7_BE29023F2D15

#include <vector>
#include <string>

extern "C" {
#include "cafe.h"
#include "family.h"
#include "cafe_shell.h"
#include "gmatrix.h"
}

class Globals;

enum LAMBDA_TYPE { UNDEFINED_LAMBDA, SINGLE_LAMBDA, MULTIPLE_LAMBDAS };

struct lambda_range
{
	double start;
	double step;
	double end;
};

struct lambda_arg_base
{
	bool search;
	LAMBDA_TYPE lambda_type;
	int bdone;
	std::string name;
	std::vector<double> lambdas;
	std::vector<double> k_weights;
	pTree lambda_tree;
	bool checkconv;
	int num_params;
	int fixcluster0;

	lambda_arg_base() : search(false), lambda_type(UNDEFINED_LAMBDA), bdone(0), 
		lambda_tree(NULL), checkconv(false), num_params(0), fixcluster0(0)
	{
	}

	virtual void load(std::vector<Argument> pargs);
	virtual const char* command() = 0;
	virtual const char* args() = 0;
	void validate_parameter_count(int expected);

	virtual int get_num_params() const
	{
		return (lambdas.size()*(k_weights.size() - fixcluster0)) + (k_weights.size() - 1);
	}
};

struct lambda_args : lambda_arg_base
{
	lambda_args() : lambda_arg_base(), vlambda(0.0), each(false), write_files(false) {}
	virtual void load(std::vector<Argument> pargs);
	double vlambda;
	std::string outfile;
	bool each;
	bool write_files;
	std::vector<lambda_range> range;

	virtual const char* command() { return "lambda"; }
	virtual const char* args() { return "lambdas (-l)"; }
};

int cafe_cmd_lambda(Globals& globals, std::vector<std::string> tokens);
void set_all_lambdas(pCafeParam param, double value);
pGMatrix cafe_lambda_distribution(pCafeParam param, const std::vector<lambda_range>& range);
void cafe_set_prior_rfsize_poisson_lambda(std::vector<double>& prior_rfsize, int shift, double* lambda);
double cafe_set_prior_rfsize_empirical(pCafeParam param, std::vector<double>& prior_rfsize);
double get_posterior(pCafeFamily pfamily, pCafeTree pcafe, family_size_range*range, std::vector<double>& ML, std::vector<double>& MAP, std::vector<double>& prior_rfsize, int quiet);
struct posterior
{
    double max_likelihood;
    double max_posterior;
};

posterior compute_posterior(pCafeFamilyItem pitem, pCafeTree pcafe, const std::vector<double>& prior_rfsize);

struct poisson_lambda
{
    double *parameters;
    int num_params;
    int num_iterations;
    double score;
};

poisson_lambda find_poisson_lambda(pCafeFamily pfamily);

const int INIT_PARAMS = 1;
const int INIT_KWEIGHTS = 2;
void initialize_params_and_k_weights(pCafeParam param, int what);
void set_parameters(pCafeParam param, lambda_args& params);
void lambda_set(pCafeParam param, lambda_args& params);
#endif

