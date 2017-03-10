#ifndef LAMBDA_H_5CC2A46B_4881_4C8F_BDE7_BE29023F2D15
#define LAMBDA_H_5CC2A46B_4881_4C8F_BDE7_BE29023F2D15

#include <vector>
#include <string>

extern "C" {
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

struct lambda_args
{
	bool search;
	LAMBDA_TYPE lambda_type;
	double vlambda;
	std::string outfile;
	int bdone;
	bool each;
	bool write_files;
	std::string name;
	std::vector<lambda_range> range;
	std::vector<double> lambdas;
	std::vector<double> k_weights;
	pTree lambda_tree;
	bool checkconv;
	int num_params;
	int fixcluster0;

	lambda_args() : search(false), lambda_type(UNDEFINED_LAMBDA), vlambda(0.0), bdone(0), each(false),
		write_files(false), lambda_tree(NULL), checkconv(false), num_params(0), fixcluster0(0)
	{
	}

};

lambda_args get_arguments(std::vector<Argument> pargs);
int cafe_cmd_lambda(Globals& globals, std::vector<std::string> tokens);
void set_all_lambdas(pCafeParam param, double value);
pGMatrix cafe_lambda_distribution(pCafeParam param, const std::vector<lambda_range>& range);

const int INIT_PARAMS = 1;
const int INIT_KWEIGHTS = 2;
void initialize_params_and_k_weights(pCafeParam param, int what);
void set_parameters(pCafeParam param, lambda_args& params);
void lambda_set(pCafeParam param, lambda_args& params);

#endif

