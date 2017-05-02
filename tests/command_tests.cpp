#include <stdexcept>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include "CppUTest/TestHarness.h"
#include "CppUTest/CommandLineTestRunner.h"
#include "cafe_commands.h"
#include "reports.h"
#include "Globals.h"
#include "viterbi.h"
#include "gene_family.h"

extern "C" {
#include "cafe.h"
};

using namespace std;

static pCafeTree create_tree()
{
	const char *newick_tree = "(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)";
	char tree[100];
	strcpy(tree, newick_tree);
	family_size_range range;
	range.min = range.root_min = 0;
	range.max = range.root_max = 15;
	return cafe_tree_new(tree, &range, 0, 0);
}

TEST_GROUP(CommandTests)
{
	vector<string> tokens;
	Globals globals;

	void setup()
	{
		srand(10);
		tokens.clear();
		globals.param.pcafe = NULL;
		globals.param.root_dist = NULL;
		globals.param.pfamily = NULL;
		globals.param.quiet = 1;
		globals.param.lambda = NULL;
		globals.param.str_log = NULL;
		globals.param.flog = stdout;
	}
};

TEST(CommandTests, command_list)
{
	std::map<string, cafe_command2> dispatcher = get_dispatcher();
	CHECK(dispatcher.find("viterbi") != dispatcher.end());

}
TEST(CommandTests, Test_cafe_cmd_source_prereqs)
{
	tokens.push_back("source");
	try
	{
		cafe_cmd_source(globals, tokens);
		FAIL("No exception was thrown");
	}
	catch (const runtime_error& err)
	{
		STRCMP_EQUAL("Usage: source <file>\n", err.what());
	}

	tokens.push_back("nonexistent");
	try
	{
		cafe_cmd_source(globals, tokens);
		FAIL("No exception was thrown");
	}
	catch (const io_error& err)
	{	
		STRCMP_EQUAL("ERROR (source): Cannot open nonexistent", err.what());
	}
};

TEST(CommandTests, cafe_cmd_generate_random_family)
{
	tokens.push_back("genfamily");
	CHECK_THROWS(std::runtime_error, cafe_cmd_generate_random_family(globals, tokens));
	tokens.push_back("filename");
	CHECK_THROWS(std::runtime_error, cafe_cmd_generate_random_family(globals, tokens));	// no tree
	globals.param.pcafe = create_tree();
	CHECK_THROWS(std::runtime_error, cafe_cmd_generate_random_family(globals, tokens));	// no family or root dist
}

TEST(CommandTests, cafe_cmd_date)
{
	char outbuf[10000];
	globals.param.flog = fmemopen(outbuf, 999, "w");
	cafe_cmd_date(globals, tokens);
	STRCMP_CONTAINS("2017", outbuf);	// this will start to fail in 2018
	fclose(globals.param.flog);
}

TEST(CommandTests, cafe_cmd_echo)
{
	char outbuf[10000];
	globals.param.flog = fmemopen(outbuf, 999, "w");
	tokens.push_back("echo");
	tokens.push_back("quick");
	tokens.push_back("brown");
	tokens.push_back("fox");
	cafe_cmd_echo(globals, tokens);
	STRCMP_EQUAL(" quick brown fox\n", outbuf);
	fclose(globals.param.flog);
}

TEST(CommandTests, cafe_cmd_exit)
{
	// all of these are values that could potentially be freed on exit
	pCafeParam param = &globals.param;
	param->str_log = NULL;
	param->lambda_tree = NULL;
	input_values_construct(&param->input, 10);
	param->pfamily = NULL;
	param->pcafe = NULL;
	param->prior_rfsize = NULL;
	param->MAP = NULL;
	param->ML = (double *)memory_new(10, sizeof(double));;
	param->str_fdata = NULL;
	globals.viterbi->viterbiPvalues = NULL;
	globals.viterbi->cutPvalues = NULL;

	cafe_cmd_exit(globals, tokens);

	LONGS_EQUAL(0, param->input.parameters);
	LONGS_EQUAL(0, param->ML);
}

TEST(CommandTests, cafe_command_report_prereqs)
{
	CHECK_THROWS(std::runtime_error, cafe_cmd_report(globals, tokens));

	CafeFamily fam;
	globals.param.pfamily = &fam;
	CHECK_THROWS(std::runtime_error, cafe_cmd_report(globals, tokens));

	CafeTree tree;
	globals.param.pcafe = &tree;
	CHECK_THROWS(std::runtime_error, cafe_cmd_report(globals, tokens));
}

void assert_gainloss_exception(Globals& globals, std::string expected)
{
	try
	{
		cafe_cmd_gainloss(globals, std::vector<std::string>());
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& e)
	{
		STRCMP_EQUAL(expected.c_str(), e.what());

	}
}

TEST(CommandTests, cafe_cmd_gainloss_exceptions)
{
	assert_gainloss_exception(globals, "ERROR: The gene families were not loaded. Please load gene families with the 'load' command.\n");

	CafeFamily fam;
	globals.param.pfamily = &fam;
	assert_gainloss_exception(globals, "ERROR: The tree was not loaded. Please load a tree with the 'tree' command.\n");

	CafeTree tree;
	globals.param.pcafe = &tree;
	assert_gainloss_exception(globals, "ERROR: Lambda values were not set. Please set lambda values with the 'lambda' or 'lambdamu' commands.\n");
}

TEST(CommandTests, cafe_cmd_log)
{
	tokens.push_back("log");
	tokens.push_back("stdout");
	cafe_cmd_log(globals, tokens);
	LONGS_EQUAL(stdout, globals.param.flog);

	tokens[1] = "log.txt";
	cafe_cmd_log(globals, tokens);
	STRCMP_EQUAL("log.txt", globals.param.str_log->buf);
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
		cafe_cmd_load(globals, tokens);
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
		cafe_cmd_load(globals, tokens);
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
		cafe_cmd_save(globals, tokens);
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

	globals.param.pcafe = NULL;
	globals.param.old_branchlength = NULL;
	cafe_cmd_tree(globals, tokens);
	CHECK(globals.param.pcafe != NULL);
	LONGS_EQUAL(8, globals.param.num_branches);
	CHECK(globals.param.old_branchlength != NULL);
	LONGS_EQUAL(212, globals.param.sum_branch_length);
	LONGS_EQUAL(81, globals.param.max_branch_length);
}

TEST(CommandTests, cafe_cmd_tree_syncs_family_if_loaded)
{
	tokens.push_back("tree");
	tokens.push_back("(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)");

	globals.param.pcafe = NULL;
	globals.param.old_branchlength = NULL;

	const char *species[] = { "", "", "chimp", "human", "mouse", "rat", "dog" };
	globals.param.pfamily = cafe_family_init(build_arraylist(species, 7));
	const char *values[] = { "description", "id", "3", "5", "7", "11", "13" };
	cafe_family_add_item(globals.param.pfamily, build_arraylist(values, 7));

	LONGS_EQUAL(-1, globals.param.pfamily->index[0]);
	cafe_cmd_tree(globals, tokens);
	LONGS_EQUAL(0, globals.param.pfamily->index[0]);
	LONGS_EQUAL(2, globals.param.pfamily->index[1]);
	LONGS_EQUAL(4, globals.param.pfamily->index[2]);
}

TEST(CommandTests, cafe_cmd_tree_missing_branch_length)
{
	tokens.push_back("tree");
	tokens.push_back("(((chimp:6,human):81,(mouse:17,rat:17):70):6,dog:9)");

	try
	{
		cafe_cmd_tree(globals, tokens);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& e)
	{
		STRCMP_EQUAL("Failed to load tree from provided string (branch length missing)", e.what());
	}
}

extern pBirthDeathCacheArray probability_cache;

void prepare_viterbi(CafeParam& param)
{
	param.pcafe = create_tree();
	double lambdas[] = { 1.5, 2.5, 3.5 };
	param.lambda = lambdas;

	const char *species[] = { "", "", "chimp", "human", "mouse", "rat", "dog" };
	param.pfamily = cafe_family_init(build_arraylist(species, 7));
	const char *values[] = { "description", "id", "3", "5", "7", "11", "13" };
	cafe_family_add_item(param.pfamily, build_arraylist(values, 7));

	cafe_family_set_species_index(param.pfamily, param.pcafe);

	param.family_size.min = param.family_size.root_min = 0;
	param.family_size.max = param.family_size.root_max = 15;

	probability_cache = birthdeath_cache_init(20);
}

TEST(CommandTests, cafe_cmd_viterbi_id_not_existing)
{
	prepare_viterbi(globals.param);

	tokens.push_back("viterbi");
	tokens.push_back("-id");
	tokens.push_back("fish");

	try
	{
		cafe_cmd_viterbi(globals, tokens);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& e)
	{
		STRCMP_EQUAL("ERROR(viterbi): fish not found", e.what());
	}
}

TEST(CommandTests, cafe_cmd_viterbi_family_out_of_range)
{
	prepare_viterbi(globals.param);

	tokens.push_back("viterbi");
	tokens.push_back("-idx");
	tokens.push_back("1000");

	try
	{
		cafe_cmd_viterbi(globals, tokens);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& e)
	{
		STRCMP_EQUAL("ERROR(viterbi): Out of range[0~1]: 1000", e.what());
	}
}

TEST(CommandTests, cafe_cmd_viterbi_args)
{
	tokens.push_back("viterbi");
	tokens.push_back("-id");
	tokens.push_back("fish");

	struct viterbi_args args = get_viterbi_arguments(build_argument_list(tokens));
	STRCMP_EQUAL("fish", args.item_id.c_str());

	tokens[1] = "-idx";
	tokens[2] = "4";
	args = get_viterbi_arguments(build_argument_list(tokens));
	LONGS_EQUAL(4, args.idx);

	tokens[1] = "-all";
	tokens[2] = "vit.txt";
	args = get_viterbi_arguments(build_argument_list(tokens));
	STRCMP_EQUAL("vit.txt", args.file.c_str());
}

TEST(CommandTests, viterbi_write)
{
	prepare_viterbi(globals.param);
	ostringstream ost;
	viterbi_write(ost, globals.param.pcafe, globals.param.pfamily);
	STRCMP_CONTAINS("id\t0\t(((chimp_3:6,human_5:6)_0:81,(mouse_7:17,rat_11:17)_0:70)_0:6,dog_13:9)_0\t0\n", ost.str().c_str());
	STRCMP_CONTAINS("Score: -inf\n", ost.str().c_str());


}

TEST(CommandTests, get_pvalue_arguments)
{
	tokens.push_back("pvalue");
	pvalue_args args = get_pvalue_arguments(build_argument_list(tokens));
	CHECK(args.infile.empty());
	CHECK(args.outfile.empty());
	LONGS_EQUAL(-1, args.index);

	tokens.push_back("-i");
	tokens.push_back("infile");
	tokens.push_back("-o");
	tokens.push_back("outfile");
	tokens.push_back("-idx");
	tokens.push_back("17");

	args = get_pvalue_arguments(build_argument_list(tokens));
	STRCMP_EQUAL("infile", args.infile.c_str());
	STRCMP_EQUAL("outfile", args.outfile.c_str());
	LONGS_EQUAL(17, args.index);
}

TEST(CommandTests, get_lhtest_arguments)
{
	tokens.push_back("lhtest");
	lhtest_args args = get_lhtest_arguments(build_argument_list(tokens));
	CHECK(args.directory.empty());
	CHECK(args.outfile.empty());
	CHECK(args.tree.empty());
	DOUBLES_EQUAL(0.0, args.lambda, .0001);

	tokens.push_back("-t");
	tokens.push_back("atree");
	tokens.push_back("-o");
	tokens.push_back("outfile");
	tokens.push_back("-l");
	tokens.push_back("0.03");
	tokens.push_back("-d");
	tokens.push_back("directory");

	args = get_lhtest_arguments(build_argument_list(tokens));
	STRCMP_EQUAL("atree", args.tree.c_str());
	STRCMP_EQUAL("outfile", args.outfile.c_str());
	STRCMP_EQUAL("directory", args.directory.c_str());
	DOUBLES_EQUAL(0.03, args.lambda, 0.00001);
}

TEST(CommandTests, cafe_cmd_lhtest)
{
	tokens.push_back("lhtest");
	try
	{
		cafe_cmd_lhtest(globals, tokens);
		FAIL("Expected exception not thrown");
	}
	catch (std::runtime_error& e)
	{
		STRCMP_EQUAL("Failed to read directory", e.what());
	}
}

TEST(CommandTests, cafe_cmd_esterror)
{
	tokens.push_back("esterror");
	esterror_args args = get_esterror_arguments(build_argument_list(tokens));
	LONGS_EQUAL(2, args.max_diff);
	CHECK_FALSE(args.peakzero);
	CHECK_FALSE(args.symmetric);
	tokens.push_back("-symm");
	tokens.push_back("-peakzero");
	tokens.push_back("-diff");
	tokens.push_back("10");
	args = get_esterror_arguments(build_argument_list(tokens));
	LONGS_EQUAL(10, args.max_diff);
	CHECK(args.peakzero);
	CHECK(args.symmetric);
}

TEST(CommandTests, validate_esterror)
{
	esterror_args args;
	CHECK_THROWS(std::runtime_error, validate(args));	// must have an output file

	args.outfile = "outfile.txt";
	CHECK_THROWS(std::runtime_error, validate(args));	// must have some data error files

	args.data_error_files.push_back("error.txt");
	CHECK_THROWS(std::runtime_error, validate(args));	// must have a truth file with one data error file

	args.truth_file = "outfile.txt";
	validate(args);	// woohoo!

	args.truth_file.clear();
	args.data_error_files.push_back("error2.txt");
	validate(args);	// or two error files  is okay too
}

