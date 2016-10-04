#ifndef LAMBDA_H_5CC2A46B_4881_4C8F_BDE7_BE29023F2D15
#define LAMBDA_H_5CC2A46B_4881_4C8F_BDE7_BE29023F2D15

#include <vector>
#include <string>
#include "family.h"

struct lambda_args
{
	bool search;
	int blambda;
	double vlambda;
	CafeParam tmp_param;
	pArgument pdist;
	pArgument pout;
	int bdone;
	bool each;
	bool write_files;
	std::string name;

	lambda_args() : search(false), blambda(0), vlambda(0.0), tmp_param(), pdist(NULL), pout(NULL), bdone(0), each(false),
		write_files(false)
	{
		memset(&tmp_param, 0, sizeof(CafeParam));
		tmp_param.posterior = 1;
	}
};

lambda_args get_arguments(pArrayList pargs);
pArrayList lambda_build_argument(std::vector<std::string> tokens);
int cafe_cmd_lambda(std::vector<std::string> tokens);

#endif

