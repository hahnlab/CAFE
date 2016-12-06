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
#include <viterbi.h>
};

#include <cafe_commands.h>
#include <lambda.h>
#include <reports.h>
#include <likelihood_ratio.h>

extern "C" {
	void show_sizes(FILE*, pCafeParam param, pCafeFamilyItem pitem, int i);
	void phylogeny_lambda_parse_func(pTree ptree, pTreeNode ptnode);
	extern pBirthDeathCacheArray probability_cache;
}


static void init_cafe_tree()
{
	const char *newick_tree = "(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)";
	char buf[100];
	strcpy(buf, "tree ");
	strcat(buf, newick_tree);
	cafe_shell_dispatch_command(buf);
}

static pCafeTree create_tree()
{
	const char *newick_tree = "(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)";
	char tree[100];
	strcpy(tree, newick_tree);
	int family_sizes[2] = { 0,15 };
	int rootfamily_sizes[2] = { 0,15 };
	return cafe_tree_new(tree, family_sizes, rootfamily_sizes, 0, 0);
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

TEST_GROUP(LikelihoodRatio)
{
	void setup()
	{
		srand(10);
	}
};

TEST_GROUP(TreeTests)
{
	void setup()
	{
		srand(10);
	}
};

TEST_GROUP(ReportTests)
{
	void setup()
	{
		srand(10);
	}
};


TEST(TreeTests, node_set_birthdeath_matrix)
{
	std::string str;

	pBirthDeathCacheArray cache = birthdeath_cache_init(10);
	pTree tree = (pTree)create_tree();
	pCafeNode node = (pCafeNode)tree_get_child(tree->root, 0);

	POINTERS_EQUAL(NULL, node->birthdeath_matrix);

	// if branch length is not set, no probabilities can be set
	node->super.branchlength = -1;
	node_set_birthdeath_matrix(node, cache, 0);
	POINTERS_EQUAL(NULL, node->birthdeath_matrix);

	// if param_lambdas not set, node's birthdeath matrix will be set
	node->super.branchlength = 6;
	node_set_birthdeath_matrix(node, cache, 0);
	CHECK(node->birthdeath_matrix != NULL)

		// if param_lambdas not set, node's birthdeath matrix will be set
	node->super.branchlength = 6;
	node_set_birthdeath_matrix(node, cache, 5);
	CHECK(node->birthdeath_matrix != NULL);

	// even if param_lambdas is set, node's birthdeath matrix will be set if num_lambdas is 0
	node->birthdeath_matrix = NULL;
	node->birth_death_probabilities.param_lambdas = (double*)memory_new(5, sizeof(double));
	node_set_birthdeath_matrix(node, cache, 0);
	CHECK(node->birthdeath_matrix != NULL);

	// if param_lambdas is set and num_lambdas > 0, put the matrices into k_bd
	node->birthdeath_matrix = NULL;
	node->k_bd = arraylist_new(5);
	node_set_birthdeath_matrix(node, cache, 5);
	POINTERS_EQUAL(NULL, node->birthdeath_matrix);
	LONGS_EQUAL(5, node->k_bd->size);
}

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

TEST(TreeTests, TestCafeTree)
{
	pCafeTree cafe_tree = create_tree();
	LONGS_EQUAL(104, cafe_tree->super.size); // tracks the size of the structure for copying purposes, etc.

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

	param.ML = (double*)memory_new(15, sizeof(double));
	param.MAP = (double*)memory_new(15, sizeof(double));

	param.family_sizes[0] = param.rootfamily_sizes[0] = 0;
	param.family_sizes[1] = param.rootfamily_sizes[1] = 15;

	LONGS_EQUAL(16,	param.pcafe->size_of_factor);	// as a side effect of create_tree

	pArrayList node_list = param.pcafe->super.nlist;
	pTreeNode node = (pTreeNode)node_list->array[1];
	CHECK(node->children->head != NULL);

	reset_birthdeath_cache(param.pcafe, 0, param.family_sizes, param.rootfamily_sizes);

	DOUBLES_EQUAL(-1.0, cafe_get_posterior(&param), 0.01);	// -1 represents an error - empirical posterior not defined. Is this safe?

	cafe_set_prior_rfsize_empirical(&param);
	CHECK_FALSE(isfinite(cafe_get_posterior(&param)));
};

TEST(TreeTests, compute_internal_node_likelihoode)
{
	pCafeTree pcafe = create_tree();
	pCafeNode node = (pCafeNode)pcafe->super.nlist->array[3];
	probability_cache = birthdeath_cache_init(pcafe->size_of_factor);
	compute_internal_node_likelihood((pTree)pcafe, (pTreeNode)node);
	DOUBLES_EQUAL(0, node->likelihoods[0], .001);
}

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
	DOUBLES_EQUAL(0.0, param.prior_rfsize[0], .001);
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

TEST(TreeTests, cafe_tree_random_probabilities)
{
	int num_families = 1;
	pCafeTree tree = create_tree();
	probability_cache = (pBirthDeathCacheArray)memory_new(num_families, sizeof(BirthDeathCacheArray));
	probability_cache->maxFamilysize = num_families;
	pArrayList node_list = tree->super.nlist;
	square_matrix bd;
	square_matrix_init(&bd, tree->familysizes[1] + 1);
	for (int i = 0; i < num_families+1; ++i)
	{
		square_matrix_set(&bd, i, 0, i);
	}
	for (int i = 0; i < node_list->size; ++i)
	{
		pCafeNode node = (pCafeNode)arraylist_get(node_list, i);
		node->birthdeath_matrix = &bd;
	}

	double *trials = cafe_tree_random_probabilities(tree, 1, 5);
	DOUBLES_EQUAL(0.0, trials[0], .001);
	DOUBLES_EQUAL(0.0, trials[1], .001);
	DOUBLES_EQUAL(0.0, trials[2], .001);
	DOUBLES_EQUAL(0.0, trials[3], .001);
	DOUBLES_EQUAL(0.0, trials[4], .001);
}

TEST(ReportTests, get_report_parameters)
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

TEST(ReportTests, write_report)
{
	std::ostringstream ost;
	CafeParam param;
	CafeFamily fam;
	fam.flist = arraylist_new(1);

	param.viterbi.num_nodes = 0;
	param.pcafe = create_tree();
	param.pfamily = &fam;
	probability_cache = NULL;
	param.num_lambdas = 3;
	double lambdas[] = { 1.5, 2.5, 3.5 };
	param.lambda = lambdas;
	param.lambda_tree = (pTree)param.pcafe;
	cafe_report(&param, ost);
	STRCMP_CONTAINS("Tree:(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)\n", ost.str().c_str());
	STRCMP_CONTAINS("Lambda:\t1.5\t2.5\t3.5\n", ost.str().c_str());
	STRCMP_CONTAINS("Lambda tree:\t(((:6,:6):81,(:17,:17):70):6,:9)\n", ost.str().c_str());
	STRCMP_CONTAINS("IDs of nodes:(((chimp<0>,human<2>)<1>,(mouse<4>,rat<6>)<5>)<3>,dog<8>)<7>\n", ost.str().c_str());
	STRCMP_CONTAINS("# Output format for: ' Average Expansion', 'Expansions', 'No Change', 'Contractions', and 'Branch-specific P-values' = (node ID, node ID): ", ost.str().c_str());
	STRCMP_CONTAINS("(0,2) (1,5) (4,6) (3,8) \n", ost.str().c_str());
	STRCMP_CONTAINS("# Output format for 'Branch cutting P-values' and 'Likelihood Ratio Test': (0, 1, 2, 3, 4, 5, 6, 7, 8)\n", ost.str().c_str());
	STRCMP_CONTAINS("", ost.str().c_str());
}

TEST(ReportTests, write_viterbi)
{
	std::ostringstream ost;
	viterbi_parameters viterbi;
	viterbi.num_nodes = 6;
	double expansion[] = { 1.3, 2.3, 3.3, 4.3, 5.3, 6.3 };
	viterbi.averageExpansion = expansion;
	int* expa[3];
	int a2[] = { 1, 3, 5, 7, 11, 13 };
	int a3[] = { 17, 19, 23, 29, 31, 37 };
	int a4[] = { 41, 43, 47, 53, 59, 61 };
	expa[0] = a2;
	expa[1] = a3;
	expa[2] = a4;
	viterbi.expandRemainDecrease = (int **)expa;
	write_viterbi(ost, viterbi);
	STRCMP_CONTAINS("Average Expansion:\t(1.3,2.3)\t(3.3,4.3)\t(5.3,6.3)\n", ost.str().c_str());
	STRCMP_CONTAINS("Expansion :\t(1,3)\t(5,7)\t(11,13)\n", ost.str().c_str());
	STRCMP_CONTAINS("Remain :\t(17,19)\t(23,29)\t(31,37)\n", ost.str().c_str());
	STRCMP_CONTAINS("Decrease :\t(41,43)\t(47,53)\t(59,61)\n", ost.str().c_str());
}


TEST(ReportTests, write_families_header)
{
	std::ostringstream ost;
	write_families_header(ost, NULL, NULL);
	STRCMP_EQUAL("'ID'\t'Newick'\t'Family-wide P-value'\t'Viterbi P-values'\n", ost.str().c_str());

	std::ostringstream ost2;
	write_families_header(ost2, (double **)1, NULL);
	STRCMP_EQUAL("'ID'\t'Newick'\t'Family-wide P-value'\t'Viterbi P-values'\t'cut P-value'\n", ost2.str().c_str());

	std::ostringstream ost3;
	write_families_header(ost3, (double **)1, (double **)1);
	STRCMP_EQUAL("'ID'\t'Newick'\t'Family-wide P-value'\t'Viterbi P-values'\t'cut P-value'\t'Likelihood Ratio'\n", ost3.str().c_str());
}

TEST(ReportTests, write_families_line)
{
	CafeParam param;
	param.likelihoodRatios = NULL;
	pCafeTree tree = create_tree();
	param.pcafe = tree;
	const char *species[] = { "", "", "chimp", "human", "mouse", "rat", "dog" };
	param.pfamily = cafe_family_init(build_arraylist(species, 7));
	cafe_family_set_species_index(param.pfamily, tree);
	const char *values[] = { "description", "id", "3", "5", "7", "11", "13" };
	cafe_family_add_item(param.pfamily, build_arraylist(values, 7));

	param.viterbi.num_nodes = 6;
	double expansion[] = { 1.3, 2.3, 3.3, 4.3, 5.3, 6.3 };
	param.viterbi.viterbiPvalues = (double**)memory_new_2dim(6, 1, sizeof(double));
	param.viterbi.viterbiPvalues[0][0] = .025;
	param.viterbi.viterbiNodeFamilysizes = (int**)memory_new_2dim(6, 1, sizeof(int));
	param.viterbi.cutPvalues = NULL;

	double maxP = .1;
	param.viterbi.maximumPvalues = &maxP;

	std::ostringstream ost;
	write_families_line(ost, &param, 0, "NodeZero");
	STRCMP_EQUAL("NodeZero\t(((chimp_3:6,human_5:6)_0:81,(mouse_7:17,rat_11:17)_0:70)_0:6,dog_13:9)_0\t0.1\t((0.025,0),(0,0),(0,0))\t\n", ost.str().c_str());

	param.viterbi.cutPvalues = (double**)memory_new_2dim(6, 1, sizeof(double));
	std::ostringstream ost2;
	write_families_line(ost2, &param, 0, "NodeZero");
	STRCMP_EQUAL("NodeZero\t(((chimp_3:6,human_5:6)_0:81,(mouse_7:17,rat_11:17)_0:70)_0:6,dog_13:9)_0\t0.1\t((0.025,0),(0,0),(0,0))\t(0,0,0,0,0,0)\t\n", ost2.str().c_str());

	param.viterbi.cutPvalues = NULL;
	param.likelihoodRatios = (double**)memory_new_2dim(tree->super.size, 1, sizeof(double));
	std::ostringstream ost3;
	write_families_line(ost3, &param, 0, "NodeZero");
	STRCMP_EQUAL("NodeZero\t(((chimp_3:6,human_5:6)_0:81,(mouse_7:17,rat_11:17)_0:70)_0:6,dog_13:9)_0\t0.1\t((0.025,0),(0,0),(0,0))\t(0,0,0,0,0,0,0,0,0)\n", ost3.str().c_str());
}

TEST(FirstTestGroup, cafe_tree_new_empty_node)
{
	pTree tree = (pTree)create_tree();
	pCafeNode node = (pCafeNode)cafe_tree_new_empty_node(tree);
	POINTERS_EQUAL(NULL, node->errormodel);
	POINTERS_EQUAL(NULL, node->birthdeath_matrix);
	POINTERS_EQUAL(NULL, node->k_bd);
	POINTERS_EQUAL(NULL, node->k_likelihoods);
	POINTERS_EQUAL(NULL, node->birth_death_probabilities.param_mus);
	POINTERS_EQUAL(NULL, node->birth_death_probabilities.param_lambdas);
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

TEST(FirstTestGroup, square_matrix_resize)
{
	square_matrix matrix;
	square_matrix_init(&matrix, 2);
	square_matrix_set(&matrix, 0, 0, 1);
	square_matrix_set(&matrix, 0, 1, 2);
	square_matrix_set(&matrix, 1, 0, 3);
	square_matrix_set(&matrix, 1, 1, 4);
	square_matrix_resize(&matrix, 3);
	LONGS_EQUAL(1, square_matrix_get(&matrix, 0, 0));
	LONGS_EQUAL(2, square_matrix_get(&matrix, 0, 1));
	LONGS_EQUAL(3, square_matrix_get(&matrix, 1, 0));
	LONGS_EQUAL(4, square_matrix_get(&matrix, 1, 1));
	square_matrix_resize(&matrix, 1);
	LONGS_EQUAL(1, square_matrix_get(&matrix, 0, 0));
}
TEST(FirstTestGroup, compute_birthdeath_rates)
{
	chooseln_cache_init(3);
	struct square_matrix* matrix = compute_birthdeath_rates(10, 0.02, 0.01, 3);
	LONGS_EQUAL(4, matrix->size);

	// matrix.values should be set to a 3x3 array of doubles
	DOUBLES_EQUAL(1, square_matrix_get(matrix, 0, 0), 0.001);
	DOUBLES_EQUAL(0, square_matrix_get(matrix, 0, 1), 0.001);
	DOUBLES_EQUAL(0, square_matrix_get(matrix, 0, 2), 0.001);
	DOUBLES_EQUAL(.086, square_matrix_get(matrix, 1, 0), 0.001);
	DOUBLES_EQUAL(.754, square_matrix_get(matrix, 1, 1), 0.001);
	DOUBLES_EQUAL(.131, square_matrix_get(matrix, 1, 2), 0.001);
	DOUBLES_EQUAL(.007, square_matrix_get(matrix, 2, 0), 0.001);
	DOUBLES_EQUAL(.131, square_matrix_get(matrix, 2, 1), 0.001);
	DOUBLES_EQUAL(.591, square_matrix_get(matrix, 2, 2), 0.001);
}

TEST(FirstTestGroup, clear_tree_viterbis)
{
	pCafeTree tree = create_tree();
	pCafeNode pcnode = (pCafeNode)tree->super.nlist->array[4];
	pcnode->familysize = 5;
	pcnode->viterbi[0] = 9;
	pcnode->viterbi[1] = 13;
	clear_tree_viterbis(tree);
	DOUBLES_EQUAL(0.0, pcnode->viterbi[0], 0.0001);
	DOUBLES_EQUAL(0.0, pcnode->viterbi[1], 0.0001);
	LONGS_EQUAL(0, pcnode->familysize);
}

int get_family_size(pTree ptree, int id)
{
	for (int i = 0; i < ptree->nlist->size; ++i)
	{
		pCafeNode node = (pCafeNode)ptree->nlist->array[i];
		if (node->super.super.id == id)
			return node->familysize;
	}
	return -1;
}
TEST(FirstTestGroup, cafe_tree_random_familysize)
{
	pCafeTree tree = create_tree();
	probability_cache = NULL;
	CafeParam param;
	param.pcafe = tree;
	param.family_sizes[0] = 1;
	param.rootfamily_sizes[0] = 1;
	param.family_sizes[1] = 10;
	param.rootfamily_sizes[1] = 3;

	reset_birthdeath_cache(param.pcafe, 0, param.family_sizes, param.rootfamily_sizes);
	pCafeNode node5 = (pCafeNode)tree->super.nlist->array[5];
	for (int i = 0; i <= 10; ++i)
		for (int j = 0; j <= 10; ++j)
			square_matrix_set(node5->birthdeath_matrix, i, j, .1);

	int max = cafe_tree_random_familysize(tree, 10);
	LONGS_EQUAL(10, max);
	LONGS_EQUAL(10, get_family_size((pTree)tree, 3));
	LONGS_EQUAL(10, get_family_size((pTree)tree, 7));
	LONGS_EQUAL(8, get_family_size((pTree)tree, 5));
}

TEST(FirstTestGroup, reset_birthdeath_cache)
{
	pCafeTree tree = create_tree();
	probability_cache = NULL;
	CafeParam param;
	param.pcafe = tree;
	param.family_sizes[0] = 0;
	param.rootfamily_sizes[0] = 0;
	param.family_sizes[1] = 10;
	param.rootfamily_sizes[1] = 10;
	reset_birthdeath_cache(param.pcafe, 0, param.family_sizes, param.rootfamily_sizes);
	CHECK_FALSE(probability_cache == NULL);
	// every node should have its birthdeath_matrix set to an entry in the cache
	// matching its branch length, lambda and mu values
	pCafeNode node = (pCafeNode)tree->super.nlist->array[3];
	struct square_matrix* expected = birthdeath_cache_get_matrix(probability_cache, node->super.branchlength, node->birth_death_probabilities.lambda, node->birth_death_probabilities.mu);
	POINTERS_EQUAL(expected, node->birthdeath_matrix);
}

TEST(FirstTestGroup, get_num_trials)
{
	std::vector<std::string> tokens;
	LONGS_EQUAL(1, get_num_trials(tokens));
	tokens.push_back("not much");
	LONGS_EQUAL(1, get_num_trials(tokens));
	tokens.push_back("-t");
	tokens.push_back("17");
	LONGS_EQUAL(17, get_num_trials(tokens));
}

TEST(FirstTestGroup, initialize_leaf_likelihoods)
{
	int rows = 5;
	int cols = 3;
	double **matrix = (double**)memory_new_2dim(rows, cols, sizeof(double));
	initialize_leaf_likelihoods_for_viterbi(matrix, rows, 3, 1, cols, NULL);
	double expected[5][3] = { { 0, 1, 0},{0,1,0},{0,1,0},{0,1,0},{0,1,0}};
	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < cols; ++j)
			DOUBLES_EQUAL(expected[i][j], matrix[i][j], .001);


	for (int i = 0; i < rows; ++i)
		expected[i][0] = 1;
	initialize_leaf_likelihoods_for_viterbi(matrix, rows, 2, -1, cols, NULL);
	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < cols; ++j)
			DOUBLES_EQUAL(expected[i][j], matrix[i][j], .001);
}

TEST(FirstTestGroup, get_clusters)
{
	double k_weights[3] = { 1,2,3 };
	get_clusters(3, 1, k_weights);
}

TEST(FirstTestGroup, write_node_headers)
{
	pCafeTree tree = create_tree();
	std::ostringstream ost1, ost2;
	write_node_headers(ost1, ost2, tree);
	STRCMP_EQUAL("DESC\tFID\tchimp\thuman\tmouse\trat\tdog\n", ost1.str().c_str());
	STRCMP_EQUAL("DESC\tFID\tchimp\t-1\thuman\t-3\tmouse\t-5\trat\t-7\tdog\n", ost2.str().c_str());
}

TEST(FirstTestGroup, write_version)
{
	std::ostringstream ost;
	write_version(ost);
	STRCMP_CONTAINS("Version: 3.2, built at", ost.str().c_str());
}

TEST(FirstTestGroup, compute_viterbis)
{
	square_matrix matrix;
	square_matrix_init(&matrix, 6);
	square_matrix_set(&matrix, 0, 0, 1);
	square_matrix_set(&matrix, 0, 1, 2);
	square_matrix_set(&matrix, 1, 0, 3);
	square_matrix_set(&matrix, 1, 1, 4);
	CafeNode node;
	node.k_bd = arraylist_new(5);
	arraylist_add(node.k_bd, &matrix);
	arraylist_add(node.k_bd, &matrix);
	node.k_likelihoods = NULL;
	reset_k_likelihoods(&node, 5, 5);
	node.k_likelihoods[0][0] = 5;
	node.k_likelihoods[0][1] = 6;
	node.k_likelihoods[0][2] = 7;
	node.k_likelihoods[0][3] = 8;
	double factors[5];
	int viterbis[10];
	node.viterbi = viterbis;
	compute_viterbis(&node, 0, factors, 0, 1, 0, 1);

	LONGS_EQUAL(1, node.viterbi[0]);
	LONGS_EQUAL(1, node.viterbi[1]);

	DOUBLES_EQUAL(12, factors[0], .001);
	DOUBLES_EQUAL(24, factors[1], .001);
}
void set_familysize_to_node_times_three(pTree ptree, pTreeNode pnode, va_list ap1)
{
	((pCafeNode)pnode)->familysize = pnode->id * 3;
}

TEST(FirstTestGroup, write_leaves)
{
	pCafeTree tree = create_tree();
	tree_traveral_infix((pTree)tree, set_familysize_to_node_times_three, NULL);
	std::ostringstream ost1, ost2;
	int id = 1234;
	int i = 42;
	write_leaves(ost1, tree, NULL, i, id, true);
	STRCMP_EQUAL("root42\t1234\t0\t6\t12\t18\t24\n", ost1.str().c_str());
	write_leaves(ost2, tree, NULL, i, id, false);
	STRCMP_EQUAL("root42\t1234\t0\t3\t6\t9\t12\t15\t18\t21\t24\n", ost2.str().c_str());

	std::ostringstream ost3, ost4;
	int k = 5;
	write_leaves(ost3, tree, &k, i, id, true);
	STRCMP_EQUAL("k5_root42\t1234\t0\t6\t12\t18\t24\n", ost3.str().c_str());

	write_leaves(ost4, tree, &k, i, id, false);
	STRCMP_EQUAL("k5_root42\t1234\t0\t3\t6\t9\t12\t15\t18\t21\t24\n", ost4.str().c_str());
}

TEST(LikelihoodRatio, cafe_likelihood_ratio_test)
{
	CafeParam param;
	param.flog = stdout;
	pCafeTree tree = create_tree();
	param.pcafe = tree;
	const char *species[] = { "", "", "chimp", "human", "mouse", "rat", "dog" };
	param.pfamily = cafe_family_init(build_arraylist(species, 7));
	param.num_threads = 1;
	cafe_likelihood_ratio_test(&param);
	DOUBLES_EQUAL(0, param.likelihoodRatios[0][0], .0001);
}

TEST(LikelihoodRatio, likelihood_ratio_report)
{
	pCafeTree tree = create_tree();
	const char *species[] = { "", "", "chimp", "human", "mouse", "rat", "dog" };
	pCafeFamily pfamily = cafe_family_init(build_arraylist(species, 7));
	const char *values[] = { "description", "family1", "3", "5", "7", "11", "13" };
	cafe_family_add_item(pfamily, build_arraylist(values, 7));
	cafe_family_set_species_index(pfamily, tree);

	std::vector<double> pvalues(2);
	pvalues[0] = 5;
	pvalues[1] = 7;

	std::vector<int> lambdas(1);
	std::vector<double*> lambda_cache(2);

	lambdas[0] = 0;
	double num = 3;
	lambda_cache[0] = &num;
	lambda_cache[1] = &num;

	char outbuf[1000];
	FILE* out = fmemopen(outbuf, 999, "w");

	likelihood_ratio_report(pfamily, tree, pvalues, lambdas, lambda_cache, out);

	fclose(out);
	STRCMP_EQUAL("family1\t(((chimp_3:6,human_5:6):81,(mouse_7:17,rat_11:17):70):6,dog_13:9)\t(0, 3.000000,0.000000)\t5\t0.025347\n", outbuf);
}

TEST(LikelihoodRatio, update_branchlength)
{
	int t = 5;
	pCafeTree tree = create_tree();
	int *old_branchlength = (int*)memory_new(tree->super.nlist->size-1, sizeof(int));
	old_branchlength[0] = 97;

	pPhylogenyNode pnode0 = (pPhylogenyNode)tree->super.nlist->array[0];
	LONGS_EQUAL(-1, pnode0->taxaid);
	pPhylogenyNode pnode1 = (pPhylogenyNode)tree->super.nlist->array[1];
	pnode1->taxaid = 1;
	update_branchlength(tree, (pTree)tree, 1.5, old_branchlength, &t);

	DOUBLES_EQUAL(6.0, pnode0->branchlength, 0.0001);	// not updated since taxaid < 0
	DOUBLES_EQUAL(688.5, pnode1->branchlength, 0.0001);	// equation applied since taxaid > 0
	LONGS_EQUAL(6, old_branchlength[0]);		// branchlengths are copied here
}

int main(int ac, char** av)
{
	return CommandLineTestRunner::RunAllTests(ac, av);
}
