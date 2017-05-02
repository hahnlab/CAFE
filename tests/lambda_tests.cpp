#include <stdexcept>
#include <sstream>

#include "CppUTest/TestHarness.h"
#include "CppUTest/CommandLineTestRunner.h"
#include "cafe_commands.h"
#include "lambda.h"
#include "lambdamu.h"
#include "Globals.h"
#include <gene_family.h>

extern "C" {
#include "cafe.h"
#include <cafe_shell.h>
	void cafe_shell_set_lambda(pCafeParam param, double* parameters);
	int __cafe_cmd_lambda_tree(pArgument parg);
};

static void init_cafe_tree(Globals& globals)
{
	const char *newick_tree = "(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)";

	char buf[100];
	strcpy(buf, "tree ");
	strcat(buf, newick_tree);
	cafe_shell_dispatch_command(globals, buf);
}

int lambda_cmd_helper(Globals& globals)
{
	std::vector<std::string> strs;
	strs.push_back("lambda");
	strs.push_back("-s");
	strs.push_back("-t");
	strs.push_back("(((2,2)1,(1,1)1)1,1)");
	return cafe_cmd_lambda(globals, strs);
}


TEST_GROUP(LambdaTests)
{
};

TEST(LambdaTests, cafe_cmd_lambda_fails_without_tree)
{
	Globals globals;
	CafeFamily fam;
	globals.param.pfamily = &fam;

	try
	{ 
		lambda_cmd_helper(globals);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& ex)
	{
		const char *expected = "ERROR: The tree was not loaded. Please load a tree with the 'tree' command.";
		STRCMP_CONTAINS(expected, ex.what());
	}
}

TEST(LambdaTests, cafe_cmd_lambda_fails_without_load)
{
	Globals globals;

	try
	{
		lambda_cmd_helper(globals);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& ex)
	{
		const char *expected = "ERROR: The gene families were not loaded. Please load gene families with the 'load' command";
		STRCMP_CONTAINS(expected, ex.what());
	}
}

TEST(LambdaTests, PrepareCafeParam)
{
	Globals globals;
	CafeFamily fam;
	CafeTree tree;
	globals.param.pfamily = &fam;
	globals.param.pcafe = &tree;
	globals.param.lambda_tree = 0;
	globals.mu_tree = NULL;
	globals.Prepare();
	POINTERS_EQUAL(0, globals.param.lambda);
	POINTERS_EQUAL(0, globals.param.mu);
	LONGS_EQUAL(-1, globals.param.num_lambdas);
	LONGS_EQUAL(-1, globals.param.num_mus);
	LONGS_EQUAL(0, globals.param.parameterized_k_value);
	POINTERS_EQUAL(cafe_shell_set_lambda, globals.param.param_set_func);
}

TEST(LambdaTests, TestCmdLambda)
{
	Globals globals;
	globals.param.quiet = 1;
	init_cafe_tree(globals);
	birthdeath_cache_init(2);
	char buf[100];
	strcpy(buf, "FAMILYDESC FAMILY Dog Chimp Human Mouse Rat");
	//cafe_shell_dispatch_command(globals, buf);

	pArrayList data = string_pchar_split(buf, ' ');
	globals.param.pfamily = cafe_family_init(data);

	strcpy(buf, "OTOPETRIN ENSF00000002390 7 7 13 18 7");
	data = string_pchar_split(buf, ' ');
	cafe_family_add_item(globals.param.pfamily, data);

	init_family_size(&globals.param.family_size, globals.param.pfamily->max_size);
	cafe_tree_set_parameters(globals.param.pcafe, &globals.param.family_size, 0);
	cafe_family_set_species_index(globals.param.pfamily, globals.param.pcafe);
	globals.param.ML = (double*)calloc(globals.param.pfamily->flist->size, sizeof(double));
	globals.param.MAP = (double*)calloc(globals.param.pfamily->flist->size, sizeof(double));

	LONGS_EQUAL(0, lambda_cmd_helper(globals));
};

TEST(LambdaTests, TestLambdaTree)
{
	Globals globals;
	init_cafe_tree(globals);
	char strs[2][100];
	strcpy(strs[0], "(((2,2)1,(1,1)1)1,1)");

	Argument arg;
	arg.argc = 1;
	char* argv[] = { strs[0], strs[1] };
	arg.argv = argv;
	__cafe_cmd_lambda_tree(&arg);
};

TEST(LambdaTests, Test_arguments)
{
	Globals globals;
	init_cafe_tree(globals);
	std::vector<std::string> strs;
	strs.push_back("lambda");
	strs.push_back("-t");
	strs.push_back("(((2,2)1,(1,1)1)1,1)");
	std::vector<Argument> pal = build_argument_list(strs);
	lambda_args args;
	args.load(pal);
	CHECK_FALSE(args.search);
	CHECK_FALSE(args.checkconv);
	CHECK_TRUE(args.lambda_tree != 0);
	LONGS_EQUAL(2, args.lambdas.size());
	DOUBLES_EQUAL(0, args.vlambda, .001);
	LONGS_EQUAL(UNDEFINED_LAMBDA, args.lambda_type);

	strs.push_back("-s");
	pal = build_argument_list(strs);
	args.load(pal);
	CHECK_TRUE(args.search);

	strs.push_back("-checkconv");
	pal = build_argument_list(strs);
	args.load(pal);
	CHECK_TRUE(args.checkconv);

	strs.push_back("-v");
	strs.push_back("14.6");
	pal = build_argument_list(strs);
	args.load(pal);
	DOUBLES_EQUAL(14.6, args.vlambda, .001);
	LONGS_EQUAL(SINGLE_LAMBDA, args.lambda_type);

	strs.push_back("-k");
	strs.push_back("19");
	pal = build_argument_list(strs);
	args.load(pal);
	LONGS_EQUAL(19, args.k_weights.size());

	strs.push_back("-f");
	pal = build_argument_list(strs);
	args.load(pal);
	LONGS_EQUAL(1, args.fixcluster0);
};


TEST(LambdaTests, Test_l_argument)
{
	Globals globals;
	init_cafe_tree(globals);
	std::vector<std::string> strs;
	strs.push_back("lambda");
	strs.push_back("-l");
	strs.push_back("15.6");
	strs.push_back("9.2");
	strs.push_back("21.8");
	std::vector<Argument> pal = build_argument_list(strs);
	lambda_args args;
	args.load(pal);
	LONGS_EQUAL(3, args.num_params);
	LONGS_EQUAL(MULTIPLE_LAMBDAS, args.lambda_type);
	DOUBLES_EQUAL(15.6, args.lambdas[0], .001);
	DOUBLES_EQUAL(9.2, args.lambdas[1], .001);
	DOUBLES_EQUAL(21.8, args.lambdas[2], .001);
};

TEST(LambdaTests, Test_p_argument)
{
	Globals globals;
	init_cafe_tree(globals);
	std::vector<std::string> strs;
	strs.push_back("lambda");
	strs.push_back("-p");
	strs.push_back("15.6");
	strs.push_back("9.2");
	strs.push_back("21.8");
	std::vector<Argument> pal = build_argument_list(strs);
	lambda_args args;
	args.load(pal);
	LONGS_EQUAL(3, args.num_params);
	DOUBLES_EQUAL(15.6, args.k_weights[0], .001);
	DOUBLES_EQUAL(9.2, args.k_weights[1], .001);
	DOUBLES_EQUAL(21.8, args.k_weights[2], .001);
};

TEST(LambdaTests, Test_r_argument)
{
	Globals globals;
	init_cafe_tree(globals);
	std::vector<std::string> strs;
	strs.push_back("lambda");
	strs.push_back("-r");
	strs.push_back("1:2:3");
	strs.push_back("-o");
	strs.push_back("test.txt");
	std::vector<Argument> pal = build_argument_list(strs);
	lambda_args args;
	args.load(pal);
	STRCMP_EQUAL("test.txt", args.outfile.c_str());
	LONGS_EQUAL(1, args.range.size());
	DOUBLES_EQUAL(1, args.range[0].start, .001);
	DOUBLES_EQUAL(2, args.range[0].step, .001);
	DOUBLES_EQUAL(3, args.range[0].end, .001);
};

TEST(LambdaTests, set_all_lambdas)
{
	CafeParam param;
	// shows that existing lambda values will be released
	param.lambda = (double *)memory_new(10, sizeof(double));
	param.num_lambdas = 15;
	set_all_lambdas(&param, 17.9);
	DOUBLES_EQUAL(17.9, param.lambda[0], 0.01);
	DOUBLES_EQUAL(17.9, param.lambda[10], 0.01);
	DOUBLES_EQUAL(17.9, param.lambda[14], 0.01);
};

TEST(LambdaTests, initialize_params_and_k_weights)
{
	CafeParam param;
	param.input.parameters = NULL;
	param.k_weights = NULL;
	param.num_params = 5;
	initialize_params_and_k_weights(&param, INIT_PARAMS);
	CHECK(param.input.parameters != NULL);
	CHECK(param.k_weights == NULL);

	param.input.parameters = NULL;
	param.parameterized_k_value = 5;
	initialize_params_and_k_weights(&param, INIT_KWEIGHTS);
	CHECK(param.input.parameters == NULL);
	CHECK(param.k_weights != NULL);
}

TEST(LambdaTests, set_parameters)
{
	lambda_args args;
	CafeParam param;

	double bar[] = { 6,7,8,9,10 };

	args.k_weights.resize(7);
	args.fixcluster0 = 3;
	args.lambdas.push_back(1);
	args.lambdas.push_back(2);
	copy(bar, bar + 5, args.k_weights.begin());
	args.num_params = 14;

	param.input.parameters = NULL;
	param.k_weights = NULL;
	param.num_lambdas = 1;

	set_parameters(&param, args);
	LONGS_EQUAL(7, param.parameterized_k_value);
	LONGS_EQUAL(14, param.num_params);
	DOUBLES_EQUAL(1.0, param.input.parameters[0], 0.001);
	DOUBLES_EQUAL(2.0, param.input.parameters[1], 0.001);
	DOUBLES_EQUAL(6.0, param.input.parameters[4], 0.001);
	DOUBLES_EQUAL(7.0, param.input.parameters[5], 0.001);
	DOUBLES_EQUAL(8.0, param.input.parameters[6], 0.001);
	DOUBLES_EQUAL(9.0, param.input.parameters[7], 0.001);
	DOUBLES_EQUAL(10.0, param.input.parameters[8], 0.001);
}

void mock_set_params(pCafeParam param, double* parameters)
{

}

TEST(LambdaTests, lambda_set)
{
	double bar[] = { 6,7,8,9,10 };

	lambda_args args;
	CafeParam param;
	args.k_weights.resize(7);
	args.fixcluster0 = 3;
	args.lambdas.push_back(1);
	args.lambdas.push_back(2);
	copy(bar, bar + 5, args.k_weights.begin());
	param.input.parameters = NULL;
	param.k_weights = NULL;
	param.num_lambdas = 1;
	param.param_set_func = mock_set_params;
	args.num_params = 14;

	lambda_set(&param, args);
}

static pCafeTree create_tree()
{
	family_size_range range;
	range.min = 0;
	range.max = 5;
	range.root_min = 0;
	range.root_max = 5;
	const char *newick_tree = "(((2,2)1,(1,1)1)1,1)";
	char tree[100];
	strcpy(tree, newick_tree);
	return cafe_tree_new(tree, &range, 0, 0);
}

TEST(LambdaTests, lambda_args__validate_parameter_count)
{
	lambda_args args;
	args.num_params = 5;
	args.lambda_tree = (pTree)create_tree();
	args.k_weights.resize(11);
	try
	{
		args.validate_parameter_count(7);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& ex)
	{
		STRCMP_CONTAINS("ERROR (lambda): The total number of parameters was not correct.\n", ex.what());
		STRCMP_CONTAINS("The total number of lambdas (-l) and proportions (-p) are 5", ex.what());
		STRCMP_CONTAINS("but 7 were expected\n", ex.what());
		STRCMP_CONTAINS("based on the tree (((2,2)1,(1,1)1)1,1) and the 11 clusters (-k).\n", ex.what());
	}

	try
	{
		args.k_weights.clear();
		args.validate_parameter_count(7);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& ex)
	{
		STRCMP_CONTAINS("ERROR (lambda): The total number of parameters was not correct.\n", ex.what());
		STRCMP_CONTAINS("The total number of lambdas (-l) are 5", ex.what());
		STRCMP_CONTAINS("but 7 were expected\n", ex.what());
		STRCMP_CONTAINS("based on the tree (((2,2)1,(1,1)1)1,1)\n", ex.what());
	}

	try
	{
		args.k_weights.resize(11);
		args.lambda_tree = NULL;
		args.validate_parameter_count(7);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& ex)
	{
		STRCMP_CONTAINS("ERROR (lambda): The total number of parameters was not correct.\n", ex.what());
		STRCMP_CONTAINS("The total number of lambdas (-l) and proportions (-p) are 5", ex.what());
		STRCMP_CONTAINS("but 7 were expected\n", ex.what());
	}

	try
	{
		args.k_weights.clear();
		args.validate_parameter_count(7);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& ex)
	{
		STRCMP_CONTAINS("ERROR (lambda): The total number of parameters was not correct.\n", ex.what());
		STRCMP_CONTAINS("The total number of lambdas (-l) are 5", ex.what());
		STRCMP_CONTAINS("but 7 were expected\n", ex.what());
	}
}

TEST(LambdaTests, lambdamu_args__validate_parameter_count)
{
	lambdamu_args args;
	args.num_params = 5;
	args.lambda_tree = (pTree)create_tree();
	args.k_weights.resize(11);

	try
	{
		args.validate_parameter_count(7);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& ex)
	{
		STRCMP_CONTAINS("ERROR (lambdamu): The total number of parameters was not correct.\n", ex.what());
		STRCMP_CONTAINS("The total number of lambdas (-l) and mus (-m) and proportions (-p) are 5", ex.what());
		STRCMP_CONTAINS("but 7 were expected\n", ex.what());
		STRCMP_CONTAINS("based on the tree (((2,2)1,(1,1)1)1,1) and the 11 clusters (-k).\n", ex.what());
	}

	try
	{
		args.k_weights.clear();
		args.validate_parameter_count(7);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& ex)
	{
		STRCMP_CONTAINS("ERROR (lambdamu): The total number of parameters was not correct.\n", ex.what());
		STRCMP_CONTAINS("The total number of lambdas (-l) and mus (-m) are 5", ex.what());
		STRCMP_CONTAINS("but 7 were expected\n", ex.what());
		STRCMP_CONTAINS("based on the tree (((2,2)1,(1,1)1)1,1)\n", ex.what());
	}

	try
	{
		args.lambda_tree = NULL;
		args.k_weights.resize(11);
		args.validate_parameter_count(7);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& ex)
	{
		STRCMP_CONTAINS("ERROR (lambdamu): The total number of parameters was not correct.\n", ex.what());
		STRCMP_CONTAINS("The total number of lambdas (-l) and mus (-m) and proportions (-p) are 5", ex.what());
		STRCMP_CONTAINS("but 7 were expected\n", ex.what());
	}

	try
	{
		args.k_weights.clear();
		args.validate_parameter_count(7);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& ex)
	{
		STRCMP_CONTAINS("ERROR (lambdamu): The total number of parameters was not correct.\n", ex.what());
		STRCMP_CONTAINS("The total number of lambdas (-l) and mus (-m) are 5", ex.what());
		STRCMP_CONTAINS("but 7 were expected\n", ex.what());
	}
}

TEST(LambdaTests, lambdamu_set)
{
	lambdamu_args args;
	args.lambdas.push_back(.1);
	args.lambdas.push_back(.2);
	args.mus.push_back(.3);
	args.mus.push_back(.4);
	args.num_params = 4;
	args.eqbg = 0;
	args.lambda_tree = (pTree)create_tree();

	CafeParam param;
	param.input.parameters = NULL;
	param.k_weights = NULL;
	param.num_lambdas = 2;
	param.param_set_func = mock_set_params;
	lambdamu_set(&param, args);
	LONGS_EQUAL(4, param.num_params);
	DOUBLES_EQUAL(.1, param.input.parameters[0], .001);
	DOUBLES_EQUAL(.2, param.input.parameters[1], .001);
	DOUBLES_EQUAL(.3, param.input.parameters[2], .001);
	DOUBLES_EQUAL(.4, param.input.parameters[3], .001);
}

TEST(LambdaTests, lambdamu_set_without_lambda_tree)
{
	lambdamu_args args;
	args.lambdas.push_back(.1);
	args.lambdas.push_back(.2);
	args.mus.push_back(.3);
	args.mus.push_back(.4);
	args.num_params = 2;
	args.eqbg = 0;
	args.lambda_tree = NULL;

	CafeParam param;
	param.input.parameters = NULL;
	param.k_weights = NULL;
	param.num_lambdas = 2;
	param.param_set_func = mock_set_params;
	lambdamu_set(&param, args);
	LONGS_EQUAL(2, param.num_params);
	DOUBLES_EQUAL(.1, param.input.parameters[0], .001);
	DOUBLES_EQUAL(.3, param.input.parameters[1], .001);
}

TEST(LambdaTests, lambdamu_set_with_k_weights)
{
	lambdamu_args args;
	args.lambdas.push_back(.1);
	args.mus.push_back(.3);
	args.k_weights.push_back(.5);
	args.k_weights.push_back(.6);

	args.num_params = 5;
	args.eqbg = 0;
	args.fixcluster0 = 0;
	args.lambda_tree = NULL;

	CafeParam param;
	param.input.parameters = NULL;
	param.k_weights = NULL;
	param.num_lambdas = 2;
	param.param_set_func = mock_set_params;
	lambdamu_set(&param, args);
	LONGS_EQUAL(5, param.num_params);
	DOUBLES_EQUAL(.1, param.input.parameters[0], .001);
	DOUBLES_EQUAL(.3, param.input.parameters[2], .001);
	DOUBLES_EQUAL(.5, param.input.parameters[4], .001);
}

TEST(LambdaTests, lambdamu_set_with_k_weights_and_tree)
{
	lambdamu_args args;
	args.lambdas.push_back(.1);
	args.lambdas.push_back(.2);
	args.mus.push_back(.3);
	args.mus.push_back(.4);
	args.k_weights.push_back(.5);
	args.k_weights.push_back(.6);

	args.num_params = 9;
	args.eqbg = 0;
	args.fixcluster0 = 0;
	args.lambda_tree = (pTree)create_tree();

	Globals globals;
	globals.param.num_lambdas = 2;
	lambdamu_set(&globals.param, args);
	LONGS_EQUAL(9, globals.param.num_params);
	DOUBLES_EQUAL(.1, globals.param.input.parameters[0], .001);
	DOUBLES_EQUAL(.2, globals.param.input.parameters[1], .001);

	DOUBLES_EQUAL(.3, globals.param.input.parameters[4], .001);
	DOUBLES_EQUAL(.4, globals.param.input.parameters[5], .001);

	DOUBLES_EQUAL(.5, globals.param.input.parameters[8], .001);
}

TEST(LambdaTests, best_lambda_mu_by_fminsearch)
{
	std::ostringstream ost;
	Globals globals;
#if 0	// TODO
	best_lambda_mu_by_fminsearch(&globals.param, 1, 1, 1, ost);
#endif
}
