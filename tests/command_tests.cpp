#include <stdexcept>
#include <vector>
#include <string>
#include <sstream>
#include "CppUTest/TestHarness.h"
#include "CppUTest/CommandLineTestRunner.h"
#include "cafe_commands.h"
#include "reports.h"

extern "C" {
#include <family.h>
#include "cafe.h"
};

using namespace std;

static pCafeTree create_tree()
{
	const char *newick_tree = "(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)";
	char tree[100];
	strcpy(tree, newick_tree);
	int family_sizes[2] = { 1,2 };
	int rootfamily_sizes[2] = { 1,2 };
	return cafe_tree_new(tree, family_sizes, rootfamily_sizes, 0, 0);
}

TEST_GROUP(CommandTests)
{
	vector<string> tokens;
	CafeParam param;

	void setup()
	{
		srand(10);
		tokens.clear();
		param.pcafe = NULL;
		param.root_dist = NULL;
		param.pfamily = NULL;
		param.quiet = 1;
		param.lambda = NULL;
		param.str_log = NULL;
		param.flog = NULL;
	}
};

typedef int(*cafe_command2)(pCafeParam cafe_param, vector<string>);

TEST(CommandTests, Test_cafe_cmd_source_prereqs)
{
	tokens.push_back("source");
	try
	{
		cafe_cmd_source(&param, tokens);
		FAIL("No exception was thrown");
	}
	catch (const runtime_error& err)
	{
		STRCMP_EQUAL("Usage: source <file>\n", err.what());
	}

	tokens.push_back("nonexistent");
	try
	{
		cafe_cmd_source(&param, tokens);
		FAIL("No exception was thrown");
	}
	catch (const runtime_error& err)
	{
		STRCMP_EQUAL("Error(source): Cannot open nonexistent\n", err.what());
	}
};

TEST(CommandTests, cafe_cmd_generate_random_family)
{
	tokens.push_back("genfamily");
	CHECK_THROWS(std::runtime_error, cafe_cmd_generate_random_family(&param, tokens));
	tokens.push_back("filename");
	CHECK_THROWS(std::runtime_error, cafe_cmd_generate_random_family(&param, tokens));	// no tree
	param.pcafe = create_tree();
	CHECK_THROWS(std::runtime_error, cafe_cmd_generate_random_family(&param, tokens));	// no family or root dist
}

TEST(CommandTests, cafe_cmd_date)
{
	char outbuf[10000];
	param.flog = fmemopen(outbuf, 999, "w");
	cafe_cmd_date(&param, tokens);
	STRCMP_CONTAINS("2016", outbuf);	// this will start to fail in 2017
	fclose(param.flog);
}

TEST(CommandTests, cafe_cmd_echo)
{
	char outbuf[10000];
	param.flog = fmemopen(outbuf, 999, "w");
	tokens.push_back("echo");
	tokens.push_back("quick");
	tokens.push_back("brown");
	tokens.push_back("fox");
	cafe_cmd_echo(&param, tokens);
	STRCMP_EQUAL(" quick brown fox\n", outbuf);
	fclose(param.flog);
}

TEST(CommandTests, cafe_cmd_exit)
{
	// all of these are values that could potentially be freed on exit
	pCafeParam param = (pCafeParam)memory_new(1, sizeof(CafeParam)); // exit frees this value from the heap
	param->str_log = NULL;
	param->mu_tree = NULL;
	param->lambda_tree = NULL;
	param->parameters = (double *)memory_new(10, sizeof(double));
	param->pfamily = NULL;
	param->pcafe = NULL;
	param->prior_rfsize_by_family = NULL;
	param->prior_rfsize = NULL;
	param->MAP = NULL;
	param->ML = (double *)memory_new(10, sizeof(double));;
	param->str_fdata = NULL;
	param->viterbi.viterbiPvalues = NULL;
	param->viterbi.cutPvalues = NULL;

	cafe_cmd_exit(param, tokens);

	LONGS_EQUAL(0, param->parameters);
	LONGS_EQUAL(0, param->ML);
}

TEST(CommandTests, cafe_command_report_prereqs)
{
	CHECK_THROWS(std::runtime_error, cafe_cmd_report(&param, tokens));

	CafeFamily fam;
	param.pfamily = &fam;
	CHECK_THROWS(std::runtime_error, cafe_cmd_report(&param, tokens));

	CafeTree tree;
	param.pcafe = &tree;
	CHECK_THROWS(std::runtime_error, cafe_cmd_report(&param, tokens));
}

void assert_gainloss_exception(CafeParam *param, std::string expected)
{
	try
	{
		cafe_cmd_gainloss(param, std::vector<std::string>());
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& e)
	{
		STRCMP_EQUAL(expected.c_str(), e.what());

	}
}

TEST(CommandTests, cafe_cmd_gainloss_exceptions)
{
	assert_gainloss_exception(&param, "ERROR: You did not load family: command 'load'\n");

	CafeFamily fam;
	param.pfamily = &fam;
	assert_gainloss_exception(&param, "ERROR: You did not specify tree: command 'tree'\n");

	CafeTree tree;
	param.pcafe = &tree;
	assert_gainloss_exception(&param, "ERROR: You did not set the parameters: command 'lambda' or 'lambdamu'\n");
}

TEST(CommandTests, cafe_cmd_log)
{
	tokens.push_back("log");
	tokens.push_back("stdout");
	cafe_cmd_log(&param, tokens);
	LONGS_EQUAL(stdout, param.flog);

	tokens[1] = "log.txt";
	cafe_cmd_log(&param, tokens);
	STRCMP_EQUAL("log.txt", param.str_log->buf);
}

TEST(CommandTests, get_load_arguments)
{
	vector<string> command = tokenize("load -t 1 -r 2 -p 0.05 -l log.txt -i fam.txt");
	struct load_args args = get_load_arguments(build_argument_list(command));
	LONGS_EQUAL(1, args.num_threads);
	LONGS_EQUAL(2, args.num_random_samples);
	DOUBLES_EQUAL(.05, args.pvalue, .000001);
	STRCMP_EQUAL("log.txt", args.log_file_name.c_str());
	STRCMP_EQUAL("fam.txt", args.family_file_name.c_str());
	CHECK(!args.filter);
}

TEST(CommandTests, cafe_cmd_load)
{
	try
	{
		tokens.push_back("load");
		cafe_cmd_load(&param, tokens);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& e)
	{
		STRCMP_EQUAL("Usage(load): load <family file>\n", e.what());
	}

	try
	{
		tokens.push_back("-t");
		tokens.push_back("5");
		cafe_cmd_load(&param, tokens);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& e)
	{
		STRCMP_EQUAL("ERROR(load): You must use -i option for input file\n", e.what());
	}
}

TEST(CommandTests, cafe_cmd_save)
{
	try
	{
		tokens.push_back("load");
		cafe_cmd_save(&param, tokens);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& e)
	{
		STRCMP_EQUAL("Usage(save): save filename", e.what());
	}

}

static pArrayList build_arraylist(const char *items[], int count)
{
	pArrayList psplit = arraylist_new(20);
	for (int i = 0; i < count; ++i)
	{
		char *str = (char*)memory_new(strlen(items[i]) + 1, sizeof(char));
		strcpy(str, items[i]);
		arraylist_add(psplit, str);
	}
	return psplit;
}

TEST(CommandTests, cafe_cmd_tree)
{
	tokens.push_back("tree");
	tokens.push_back("(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)");

	param.pcafe = NULL;
	param.old_branchlength = NULL;
	cafe_cmd_tree(&param, tokens);
	CHECK(param.pcafe != NULL);
	LONGS_EQUAL(8, param.num_branches);
	CHECK(param.old_branchlength != NULL);
	LONGS_EQUAL(212, param.sum_branch_length);
	LONGS_EQUAL(81, param.max_branch_length);
}

TEST(CommandTests, cafe_cmd_tree_syncs_family_if_loaded)
{
	tokens.push_back("tree");
	tokens.push_back("(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)");

	param.pcafe = NULL;
	param.old_branchlength = NULL;

	const char *species[] = { "", "", "chimp", "human", "mouse", "rat", "dog" };
	param.pfamily = cafe_family_init(build_arraylist(species, 7));
	const char *values[] = { "description", "id", "3", "5", "7", "11", "13" };
	cafe_family_add_item(param.pfamily, build_arraylist(values, 7));

	LONGS_EQUAL(-1, param.pfamily->index[0]);
	cafe_cmd_tree(&param, tokens);
	LONGS_EQUAL(0, param.pfamily->index[0]);
	LONGS_EQUAL(2, param.pfamily->index[1]);
	LONGS_EQUAL(4, param.pfamily->index[2]);
}

TEST(CommandTests, cafe_cmd_tree_missing_branch_length)
{
	tokens.push_back("tree");
	tokens.push_back("(((chimp:6,human):81,(mouse:17,rat:17):70):6,dog:9)");

	try
	{
		cafe_cmd_tree(&param, tokens);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& e)
	{
		STRCMP_EQUAL("Failed to load tree from provided string (branch length missing)", e.what());
	}
}