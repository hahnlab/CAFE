#include <stdexcept>
#include <vector>
#include <string>
#include <sstream>

#include "CppUTest/TestHarness.h"
#include "CppUTest/CommandLineTestRunner.h"
#include <math.h>

extern "C" {
#include <utils_string.h>
#include <cafe_shell.h>
#include <tree.h>
#include <cafe.h>
#include <chooseln_cache.h>
};

#include <cafe_commands.h>
#include <lambda.h>
#include <reports.h>

extern "C" {
	void show_sizes(FILE*, pCafeParam param, pCafeFamilyItem pitem, int i);
	void phylogeny_lambda_parse_func(pTree ptree, pTreeNode ptnode);
}


static void init_cafe_tree()
{
	const char *newick_tree = "(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)";
	char buf[100];
	strcpy(buf, "tree ");
	strcat(buf, newick_tree);
	cafe_shell_dispatch_command(buf);
}


int cafe_cmd_atoi(int argc, char* argv[])
{
	return atoi(argv[1]);
}

CafeShellCommand cafe_cmd_test[] =
{
	{ "atoi", cafe_cmd_atoi },
	{ NULL, NULL }
};


TEST_GROUP(FirstTestGroup)
{
	void setup()
	{
		srand(10);
	}
};

TEST(FirstTestGroup, TestStringSplitter)
{
	char c[10];
	c[0] = 0;
	LONGS_EQUAL(0, string_pchar_space_split(c)->size);

	pArrayList pArray;
	strcpy(c, "a b");
	pArray = string_pchar_space_split(c);
	LONGS_EQUAL(2, pArray->size);
	STRCMP_EQUAL("a", (char *)(pArray->array[0]));
	STRCMP_EQUAL("b", (char *)(pArray->array[1]));
}

TEST(FirstTestGroup, Tokenize)
{
	char c[10];
	c[0] = 0;
	LONGS_EQUAL(0, tokenize(c).size());
	strcpy(c, " ");
	LONGS_EQUAL(0, tokenize(c).size());

	strcpy(c, "a b\r\n");
	std::vector<std::string> arr = tokenize(c);
	LONGS_EQUAL(2, arr.size());
	STRCMP_EQUAL("a", arr[0].c_str());
	STRCMP_EQUAL("b", arr[1].c_str());
}

static pCafeTree create_tree()
{
	const char *newick_tree = "(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)";
	char tree[100];
	strcpy(tree, newick_tree);
	int family_sizes[2] = { 1,2 };
	int rootfamily_sizes[2] = { 1,2 };
	return cafe_tree_new(tree, family_sizes, rootfamily_sizes, 0, 0);
}

TEST(FirstTestGroup, TestCafeTree)
{
	pCafeTree cafe_tree = create_tree();
	LONGS_EQUAL(128, cafe_tree->super.size);

	// Find chimp in the tree after two branches of length 6,81,6
	pTreeNode ptnode = (pTreeNode)cafe_tree->super.root;
	CHECK(tree_is_root(&cafe_tree->super, ptnode));
	ptnode = (pTreeNode)tree_get_child(ptnode, 0);
	pPhylogenyNode pnode = (pPhylogenyNode)ptnode;
	LONGS_EQUAL(6, pnode->branchlength);	// root to first branch = 6

	ptnode = (pTreeNode)tree_get_child(ptnode, 0);
	pnode = (pPhylogenyNode)ptnode;
	LONGS_EQUAL(81, pnode->branchlength); // 1st to 2nd = 81

	ptnode = (pTreeNode)tree_get_child(ptnode, 0);
	pnode = (pPhylogenyNode)ptnode;
	STRCMP_EQUAL("chimp", pnode->name);
	LONGS_EQUAL(6, pnode->branchlength);  // 2nd to chimp leaf = 6
	CHECK(tree_is_leaf(ptnode));


}

TEST(FirstTestGroup, TestShellDispatcher)
{
	char c[10];
	cafe_shell_init(1);

	CafeShellCommand *old = cafe_cmd;
	cafe_cmd[0] = cafe_cmd_test[0];
	cafe_cmd[1] = cafe_cmd_test[1];

	strcpy(c, "atoi 9528");
	LONGS_EQUAL(9528, cafe_shell_dispatch_command(c));

	strcpy(c, "# a comment");
	LONGS_EQUAL(0, cafe_shell_dispatch_command(c));

	strcpy(c, "unknown");
	LONGS_EQUAL(CAFE_SHELL_NO_COMMAND, cafe_shell_dispatch_command(c));
}

TEST(FirstTestGroup, TestShowSizes)
{
	char outbuf[10000];
	setbuf(stdout, outbuf);

	CafeParam param;
	param.rootfamily_sizes[0] = 29;
	param.rootfamily_sizes[1] = 31;
	param.family_sizes[0] = 37;
	param.family_sizes[1] = 41;

	CafeFamilyItem item;
	item.ref = 14;
	CafeTree tree;
	tree.rootfamilysizes[0] = 11;
	tree.rootfamilysizes[1] = 13;
	tree.familysizes[0] = 23;
	tree.familysizes[1] = 19;
	tree.rfsize = 17;
	param.pcafe = &tree;
	FILE* in = fmemopen(outbuf, 999, "w");
	show_sizes(in, &param, &item, 7);
	fclose(in);
	STRCMP_CONTAINS(">> 7 14", outbuf);
	STRCMP_CONTAINS("Root size: 11 ~ 13 , 17", outbuf);
	STRCMP_CONTAINS("Family size: 23 ~ 19", outbuf);
	STRCMP_CONTAINS("Root size: 29 ~ 31", outbuf);
	STRCMP_CONTAINS("Family size: 37 ~ 41", outbuf);
}

TEST(FirstTestGroup, TestPhylogenyLoadFromString)
{
	char outbuf[10000];
	strcpy(outbuf, "(((1,1)1,(2,2)2)2,2)");
	cafe_shell_init(1);
	init_cafe_tree();
	pTree tree = phylogeny_load_from_string(outbuf, tree_new, phylogeny_new_empty_node, phylogeny_lambda_parse_func, 0);
	CHECK(tree != 0);
	LONGS_EQUAL(56, tree->size);

};

TEST(FirstTestGroup, Test_cafe_cmd_source)
{
	CafeParam param;

	std::vector<std::string> strs;
	strs.push_back("source");
	LONGS_EQUAL( -1, cafe_cmd_source(&param, strs));

	strs.push_back("nonexistent");
	LONGS_EQUAL(-1, cafe_cmd_source(&param, strs));
};

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

TEST(FirstTestGroup, Test_cafe_get_posterior)
{
	CafeParam param;
	param.flog = stdout;
	param.quiet = 1;
	param.prior_rfsize_by_family = NULL;
	param.prior_rfsize = NULL;
	param.pcafe = create_tree();
	const char *species[] = { "", "", "chimp", "human", "mouse", "rat", "dog" };
	param.pfamily = cafe_family_init(build_arraylist(species, 7));
	cafe_family_set_species_index(param.pfamily, param.pcafe);
	const char *values[] = { "description", "id", "3", "5", "7", "11", "13" };
	cafe_family_add_item(param.pfamily, build_arraylist(values, 7));

	param.ML = (double*)memory_new(5, sizeof(double));
	param.MAP = (double*)memory_new(5, sizeof(double));

	DOUBLES_EQUAL(-1.0, cafe_get_posterior(&param), 0.01);	// -1 represents an error - empirical posterior not defined. Is this safe?

	cafe_set_prior_rfsize_empirical(&param);
	CHECK_FALSE(isfinite(cafe_get_posterior(&param)));
};

TEST(FirstTestGroup, cafe_set_prior_rfsize_empirical)
{
	CafeParam param;
	param.flog = stdout;
	param.prior_rfsize = NULL;
	const char *species[] = { "", "", "chimp", "human", "mouse", "rat", "dog" };
	param.pfamily = cafe_family_init(build_arraylist(species, 7));
	const char *values[] = { "description", "id", "3", "5", "7", "11", "13" };
	cafe_family_add_item(param.pfamily, build_arraylist(values, 7));

	param.pcafe = create_tree();
	cafe_set_prior_rfsize_empirical(&param);
	DOUBLES_EQUAL(0.5679, param.prior_rfsize[0], .001);
}

TEST(FirstTestGroup, list_commands)
{
	std::ostringstream ost;
	list_commands(ost);
	STRCMP_CONTAINS("lambda", ost.str().c_str());
	STRCMP_CONTAINS("tree", ost.str().c_str());
	STRCMP_CONTAINS("load", ost.str().c_str());
	STRCMP_CONTAINS("branchlength", ost.str().c_str());
}

TEST(FirstTestGroup, cafe_cmd_date)
{
	CafeParam param;
	param.quiet = 1;
	char outbuf[10000];
	param.flog = fmemopen(outbuf, 999, "w");
	cafe_cmd_date(&param, std::vector<std::string>());
	STRCMP_CONTAINS("2016", outbuf);	// this will start to fail in 2017
	fclose(param.flog);
}

TEST(FirstTestGroup, cafe_cmd_echo)
{
	CafeParam param;
	param.quiet = 1;
	char outbuf[10000];
	param.flog = fmemopen(outbuf, 999, "w");
	std::vector<std::string> tokens;
	tokens.push_back("echo");
	tokens.push_back("quick");
	tokens.push_back("brown");
	tokens.push_back("fox");
	cafe_cmd_echo(&param, tokens);
	STRCMP_EQUAL(" quick brown fox\n", outbuf);
	fclose(param.flog);
}

TEST(FirstTestGroup, cafe_cmd_exit)
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
	std::vector<std::string> tokens;
	cafe_cmd_exit(param, tokens);
	LONGS_EQUAL(0, param->parameters);
	LONGS_EQUAL(0, param->ML);
}

void store_node_id(pTree ptree, pTreeNode ptnode, va_list ap1)
{
	va_list ap;
	va_copy(ap, ap1);
	std::vector<int> * ids = va_arg(ap, std::vector<int> *);
	ids->push_back(ptnode->id);
	va_end(ap);
}

TEST(FirstTestGroup, tree_traversal_prefix)
{
	std::vector<int> ids;
	pCafeTree tree = create_tree();
	tree_traveral_prefix(&tree->super, store_node_id, &ids);
	LONGS_EQUAL(7, ids[0]);
	LONGS_EQUAL(3, ids[1]);
	LONGS_EQUAL(1, ids[2]);
	LONGS_EQUAL(0, ids[3]);
	LONGS_EQUAL(2, ids[4]);
	LONGS_EQUAL(5, ids[5]);
	LONGS_EQUAL(4, ids[6]);
	LONGS_EQUAL(6, ids[7]);
	LONGS_EQUAL(8, ids[8]);
}

TEST(FirstTestGroup, tree_traversal_infix)
{
	std::vector<int> ids;
	pCafeTree tree = create_tree();
	tree_traveral_infix(&tree->super, store_node_id, &ids);
	LONGS_EQUAL(0, ids[0]);
	LONGS_EQUAL(1, ids[1]);
	LONGS_EQUAL(2, ids[2]);
	LONGS_EQUAL(3, ids[3]);
	LONGS_EQUAL(4, ids[4]);
	LONGS_EQUAL(5, ids[5]);
	LONGS_EQUAL(6, ids[6]);
	LONGS_EQUAL(7, ids[7]);
	LONGS_EQUAL(8, ids[8]);
}

TEST(FirstTestGroup, tree_traversal_postfix)
{
	std::vector<int> ids;
	pCafeTree tree = create_tree();
	tree_traveral_postfix(&tree->super, store_node_id, &ids);
	LONGS_EQUAL(0, ids[0]);
	LONGS_EQUAL(2, ids[1]);
	LONGS_EQUAL(1, ids[2]);
	LONGS_EQUAL(4, ids[3]);
	LONGS_EQUAL(6, ids[4]);
	LONGS_EQUAL(5, ids[5]);
	LONGS_EQUAL(3, ids[6]);
	LONGS_EQUAL(8, ids[7]);
	LONGS_EQUAL(7, ids[8]);
}

TEST(FirstTestGroup, cafe_tree_random_probabilities)
{
	int num_families = 1;
	pCafeTree tree = create_tree();
	tree->pbdc_array = (pBirthDeathCacheArray)memory_new(num_families, sizeof(BirthDeathCacheArray));
	tree->pbdc_array->maxFamilysize = num_families;
	pArrayList node_list = tree->super.nlist;
	double **bd = (double **)memory_new_2dim(tree->super.size, num_families, sizeof(double));
	for (int i = 0; i < tree->super.size; ++i)
	{
		bd[i][0] = i;
	}
	for (int i = 0; i < node_list->size; ++i)
	{
		pCafeNode node = (pCafeNode)arraylist_get(node_list, i);
		node->bd = bd;
	}

	double *trials = cafe_tree_random_probabilities(tree, 1, 5);
	DOUBLES_EQUAL(0.0, trials[0], .001);
	DOUBLES_EQUAL(0.0, trials[1], .001);
	DOUBLES_EQUAL(0.0, trials[2], .001);
	DOUBLES_EQUAL(0.0, trials[3], .001);
	DOUBLES_EQUAL(0.0, trials[4], .001);
}

TEST(FirstTestGroup, get_report_parameters)
{
	report_parameters params;

	std::vector<std::string> tokens;
	tokens.push_back("report");
	tokens.push_back("myreport");
	tokens.push_back("branchcutting");

	get_report_parameters(params, tokens);
	STRCMP_EQUAL("myreport", params.name.c_str());
	LONGS_EQUAL(1, params.bc);
	LONGS_EQUAL(0, params.lh);

	tokens[2] = "lh2";
	get_report_parameters(params, tokens);
	LONGS_EQUAL(1, params.lh2);
	LONGS_EQUAL(0, params.lh);
}


TEST(FirstTestGroup, cafe_command_report_prereqs)
{
	CafeParam param;
	CafeFamily fam;
	CafeTree tree;
	std::vector<std::string> tokens;
	param.pfamily = NULL;
	CHECK_THROWS(std::runtime_error, cafe_cmd_report(&param, tokens));
	param.pfamily = &fam;
	param.pcafe = NULL;
	CHECK_THROWS(std::runtime_error, cafe_cmd_report(&param, tokens));
	param.pcafe = &tree;
	param.lambda = NULL;
	CHECK_THROWS(std::runtime_error, cafe_cmd_report(&param, tokens));
}

TEST(FirstTestGroup, compute_likelihood)
{
	chooseln_cache cache = { 0,0 };
	pCafeTree tree = create_tree();
	int maxFamilySize = MAX(tree->rootfamilysizes[1], tree->familysizes[1]);
	chooseln_cache_init2(&cache, maxFamilySize);
	pCafeNode node = (pCafeNode)tree->super.root;
	compute_likelihood((pTree)tree, (pTreeNode)node, &cache);
	DOUBLES_EQUAL(0, node->likelihoods[0], 0.01);
}

TEST(FirstTestGroup, compute_internal_node_likelihood)
{
	chooseln_cache_init(50);
	pTree tree = (pTree)create_tree();
	pCafeNode node = (pCafeNode)tree->root;
	node->lambda = 0.01;
	node->mu = -1;
	int factors = ((pCafeTree)tree)->size_of_factor;
	for (int i = 0; i < factors; ++i)
	{
		((pCafeNode)tree_get_child(tree->root, 0))->likelihoods[i] = .5;
		((pCafeNode)tree_get_child(tree->root, 1))->likelihoods[i] = .5;
	}
	tree_get_child(tree->root, 1);
	compute_internal_node_likelihood(tree, tree->root, NULL);
	DOUBLES_EQUAL(0.214, node->likelihoods[0], 0.001);
	DOUBLES_EQUAL(0.192, node->likelihoods[1], 0.001);
}

TEST(FirstTestGroup, compute_leaf_node_likelihood)
{
	pTree tree = (pTree)create_tree();
	pCafeNode leaf = (pCafeNode)tree_get_child(tree->root, 0);

	// no family size set
	compute_leaf_node_likelihood(tree, (pTreeNode)leaf);
	DOUBLES_EQUAL(1.0, leaf->likelihoods[0], 0.001);
	DOUBLES_EQUAL(1.0, leaf->likelihoods[1], 0.001);

	leaf->familysize = 1;
	compute_leaf_node_likelihood(tree, (pTreeNode)leaf);
	DOUBLES_EQUAL(0.0, leaf->likelihoods[0], 0.001);
	DOUBLES_EQUAL(1.0, leaf->likelihoods[1], 0.001);
}

TEST(FirstTestGroup, cafe_tree_new_empty_node)
{
	pTree tree = (pTree)create_tree();
	pCafeNode node = (pCafeNode)cafe_tree_new_empty_node(tree);
	POINTERS_EQUAL(NULL, node->errormodel);
	POINTERS_EQUAL(NULL, node->bd);
	POINTERS_EQUAL(NULL, node->k_bd);
	POINTERS_EQUAL(NULL, node->k_likelihoods);
	POINTERS_EQUAL(NULL, node->param_mus);
	POINTERS_EQUAL(NULL, node->param_lambdas);
	LONGS_EQUAL(-1, node->familysize);
}

TEST(FirstTestGroup, chooseln_cache)
{
	struct chooseln_cache cache;
	cache.values = 0;
	CHECK_FALSE(chooseln_is_init2(&cache));
	chooseln_cache_init2(&cache, 10);
	CHECK_TRUE(chooseln_is_init2(&cache));
	LONGS_EQUAL(10, get_chooseln_cache_size2(&cache));
	DOUBLES_EQUAL(4.025, chooseln_get2(&cache, 8, 5), .001);
	DOUBLES_EQUAL(1.098, chooseln_get2(&cache, 3, 2), .001);
	DOUBLES_EQUAL(1.791, chooseln_get2(&cache, 6, 5), .001);
	DOUBLES_EQUAL(4.43, chooseln_get2(&cache, 9, 3), .001);
	chooseln_cache_free2(&cache);
}

TEST(FirstTestGroup, birthdeath_rate_with_log_alpha	)
{
	struct chooseln_cache cache;
	chooseln_cache_init2(&cache, 50);
	DOUBLES_EQUAL(0.107, birthdeath_rate_with_log_alpha(40, 42, -1.37, 0.5, &cache), .001);
	DOUBLES_EQUAL(0.006, birthdeath_rate_with_log_alpha(41, 34, -1.262, 0.4, &cache), .001);


}

TEST(FirstTestGroup, birthdeath_likelihood_with_s_c)
{
	struct chooseln_cache cache;
	chooseln_cache_init2(&cache, 50);
	DOUBLES_EQUAL(0.083, birthdeath_likelihood_with_s_c(40, 42, 0.42, 0.5, -1, &cache), .001);
	DOUBLES_EQUAL(0.023, birthdeath_likelihood_with_s_c(41, 34, 0.54, 0.4, -1, &cache), .001);
}





int main(int ac, char** av)
{
	return CommandLineTestRunner::RunAllTests(ac, av);
}
