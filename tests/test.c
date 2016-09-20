#include "CppUTest/TestHarness.h"
#include "CppUTest/CommandLineTestRunner.h"
extern "C"
{
#include <utils_string.h>
#include <cafe_shell.h>

}

const char *newick_tree = "(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)";

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
};

TEST(FirstTestGroup, TestStringSplitter)
{
	char c[10];
	pArrayList pArray;
	strcpy(c, "a b");
	pArray = string_pchar_space_split(c);
	LONGS_EQUAL(2, pArray->size);
	STRCMP_EQUAL("a", (char *)(pArray->array[0]));
	STRCMP_EQUAL("b", (char *)(pArray->array[1]));
}

TEST(FirstTestGroup, TestCafeFamilyNew)
{
	pCafeFamily fam;
	char fname[100];
	strcpy(fname, "Nonexistent.tab");
	POINTERS_EQUAL(NULL, cafe_family_new(fname, 1));

	strcpy(fname, "../example/example_data.tab");
	fam = cafe_family_new(fname, 1);
	LONGS_EQUAL(5, fam->num_species);
	STRCMP_EQUAL("Dog", fam->species[0]);
	LONGS_EQUAL(59, fam->flist->size);
	STRCMP_EQUAL("Rat", fam->species[4]);
	LONGS_EQUAL(59, fam->flist->size);
}

TEST(FirstTestGroup, TestCafeTree)
{
	char tree[100];
	strcpy(tree, newick_tree);
	int family_sizes[2] = { 1,1 };
	int rootfamily_sizes[2] = { 1,1 };
	pCafeTree cafe_tree = cafe_tree_new(tree, family_sizes, rootfamily_sizes, 0, 0);
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

}

TEST(FirstTestGroup, TestCmdLambda_FailsWithoutTree)
{
	char c[10];
	cafe_shell_init(1);
	char strs[4][30];
	strcpy(strs[0], "lambda");
	strcpy(strs[1], "-s");
	strcpy(strs[2], "-t");
	strcpy(strs[3], "(((2,2)1,(1,1)1)1,1)");
	char *argv[] = { strs[0], strs[1], strs[2], strs[3] };
	char buf[1000];
	setbuf(stderr, buf);

	LONGS_EQUAL(-1, cafe_cmd_lambda(4, argv));
	const char *expected = "ERROR(lambda): You did not specify tree: command 'tree'";
	STRCMP_CONTAINS(expected, buf);
	setbuf(stderr, NULL);
}

TEST(FirstTestGroup, TestCmdLambdaFailsWithoutLoad)
{
	cafe_shell_init(1);
	char strs[4][30];
	strcpy(strs[0], "lambda");
	strcpy(strs[1], "-s");
	strcpy(strs[2], "-t");
	strcpy(strs[3], "(((2,2)1,(1,1)1)1,1)");
	char *argv[] = { strs[0], strs[1], strs[2], strs[3] };
	char buf[1000];

	strcpy(buf, "tree ");
	strcat(buf, newick_tree);
	cafe_shell_dispatch_command(buf);

	char errbuf[1000];
	setbuf(stderr, errbuf);
	char outbuf[10000];
	setbuf(stdout, outbuf);
	LONGS_EQUAL(-1, cafe_cmd_lambda(4, argv));
	const char *expected = "ERROR(lambda): Please load family (\"load\") and cafe tree (\"tree\") before running \"lambda\" command.";
	STRCMP_CONTAINS(expected, errbuf);
	setbuf(stderr, NULL);
	setbuf(stdout, NULL);
}

TEST(FirstTestGroup, TestCmdLambda)
{
	cafe_shell_init(1);
	char strs[4][30];
	strcpy(strs[0], "lambda");
	strcpy(strs[1], "-s");
	strcpy(strs[2], "-t");
	strcpy(strs[3], "(((2,2)1,(1,1)1)1,1)");
	char *argv[] = { strs[0], strs[1], strs[2], strs[3] };

	char buf[100];
	strcpy(buf, "tree ");
	strcat(buf, newick_tree);
	cafe_shell_dispatch_command(buf);

	strcpy(buf, "load -i ../example/example_data.tab");
	cafe_shell_dispatch_command(buf);

	LONGS_EQUAL(0, cafe_cmd_lambda(4, argv));
}

extern "C" {
	void show_sizes(FILE*, pCafeParam param, pCafeFamilyItem pitem, int i);
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

int main(int ac, char** av)
{
	return CommandLineTestRunner::RunAllTests(ac, av);
}
