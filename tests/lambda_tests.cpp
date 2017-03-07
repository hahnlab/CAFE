#include <stdexcept>

#include "CppUTest/TestHarness.h"
#include "CppUTest/CommandLineTestRunner.h"
#include "cafe_commands.h"
#include "lambda.h"
#include "Globals.h"

extern "C" {
#include <cafe_shell.h>
	int __cafe_cmd_lambda_tree(pArgument parg);
	void cafe_shell_set_lambda(pCafeParam param, double* parameters);
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

TEST(LambdaTests, TestCmdLambda_FailsWithoutTree)
{
	Globals globals;

	try
	{ 
		lambda_cmd_helper(globals);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& ex)
	{
		const char *expected = "ERROR(lambda): You did not specify tree: command 'tree'";
		STRCMP_CONTAINS(expected, ex.what());
	}
}

TEST(LambdaTests, PrepareCafeParamFailsWithoutLoad)
{
	Globals globals;
	const char *newick_tree = "(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)";

	char buf[100];
	strcpy(buf, "tree ");
	strcat(buf, newick_tree);
	cafe_shell_dispatch_command(globals, buf);
	try
	{
		globals.Prepare();
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& err)
	{
		const char *expected = "ERROR(lambda): Please load family (\"load\") and cafe tree (\"tree\") before running \"lambda\" command.";
		STRCMP_EQUAL(expected, err.what());
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
	strcpy(buf, "load -i ../example/example_data.tab");
	cafe_shell_dispatch_command(globals, buf);

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
	lambda_args args = get_arguments(pal);
	CHECK_FALSE(args.search);
	CHECK_FALSE(args.tmp_param.checkconv);
	CHECK_TRUE(args.tmp_param.lambda_tree != 0);
	LONGS_EQUAL(2, args.tmp_param.num_lambdas);
	DOUBLES_EQUAL(0, args.vlambda, .001);
	LONGS_EQUAL(UNDEFINED_LAMBDA, args.lambda_type);

	strs.push_back("-s");
	pal = build_argument_list(strs);
	args = get_arguments(pal);
	CHECK_TRUE(args.search);

	strs.push_back("-checkconv");
	pal = build_argument_list(strs);
	args = get_arguments(pal);
	CHECK_TRUE(args.tmp_param.checkconv);

	strs.push_back("-v");
	strs.push_back("14.6");
	pal = build_argument_list(strs);
	args = get_arguments(pal);
	DOUBLES_EQUAL(14.6, args.vlambda, .001);
	LONGS_EQUAL(SINGLE_LAMBDA, args.lambda_type);

	strs.push_back("-k");
	strs.push_back("19");
	pal = build_argument_list(strs);
	args = get_arguments(pal);
	LONGS_EQUAL(19, args.tmp_param.parameterized_k_value);

	strs.push_back("-f");
	pal = build_argument_list(strs);
	args = get_arguments(pal);
	LONGS_EQUAL(1, args.tmp_param.fixcluster0);
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
	lambda_args args = get_arguments(pal);
	LONGS_EQUAL(3, args.tmp_param.num_params);
	LONGS_EQUAL(MULTIPLE_LAMBDAS, args.lambda_type);
	DOUBLES_EQUAL(15.6, args.tmp_param.lambda[0], .001);
	DOUBLES_EQUAL(9.2, args.tmp_param.lambda[1], .001);
	DOUBLES_EQUAL(21.8, args.tmp_param.lambda[2], .001);
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
	lambda_args args = get_arguments(pal);
	LONGS_EQUAL(3, args.tmp_param.num_params);
	DOUBLES_EQUAL(15.6, args.tmp_param.k_weights[0], .001);
	DOUBLES_EQUAL(9.2, args.tmp_param.k_weights[1], .001);
	DOUBLES_EQUAL(21.8, args.tmp_param.k_weights[2], .001);
};

TEST(LambdaTests, Test_r_argument)
{
	Globals globals;
	init_cafe_tree(globals);
	std::vector<std::string> strs;
	strs.push_back("lambda");
	strs.push_back("-r");
	strs.push_back("-o");
	strs.push_back("test.txt");
	std::vector<Argument> pal = build_argument_list(strs);
	lambda_args args = get_arguments(pal);
	STRCMP_EQUAL("-r", args.dist.opt);
	LONGS_EQUAL(1, args.out.argc);
	STRCMP_EQUAL("test.txt", args.out.argv[0]);
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

