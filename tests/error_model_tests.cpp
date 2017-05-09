#include <sstream>
#include <vector>

#include "CppUTest/TestHarness.h"

#include <error_model.h>
#include <gene_family.h>

extern "C" {
//#include <cafe_shell.h>
//#include <tree.h>
#include <cafe.h>
};

TEST_GROUP(ErrorModel)
{
	void setup()
	{
		srand(10);
	}
};

static pCafeTree create_tree()
{
	family_size_range range;

	range.min = range.root_min = 0;
	range.max = range.root_max = 15;

	const char *newick_tree = "(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)";
	char tree[100];
	strcpy(tree, newick_tree);
	return cafe_tree_new(tree, &range, 0, 0);
}

TEST(ErrorModel, get_error_model)
{
	std::vector<double> v(5, 0.2);

	CafeFamily fam;
	fam.errors = NULL;
	POINTERS_EQUAL(NULL, get_error_model(&fam, "nonexistent"));

	fam.errors = arraylist_new(10);
	POINTERS_EQUAL(NULL, get_error_model(&fam, "nonexistent"));

	char buf[10];
	ErrorStruct err;
	strcpy(buf, "SPECIES1");
	err.errorfilename = buf;
	arraylist_add(fam.errors, &err);
	POINTERS_EQUAL(NULL, get_error_model(&fam, "nonexistent"));

	POINTERS_EQUAL(&err, get_error_model(&fam, "species1"));
}

TEST(ErrorModel, error_model_write)
{
	std::ostringstream ost;
	ErrorStruct err;
	err.maxfamilysize = 4;
	err.fromdiff = 1;
	err.todiff = 5;
	err.errormatrix = (double**)memory_new_2dim(10, 10, sizeof(double));
	ost << err;

	STRCMP_CONTAINS("maxcnt:4\n", ost.str().c_str());
	STRCMP_CONTAINS("cntdiff 1 2 3 4 5\n", ost.str().c_str());
	STRCMP_CONTAINS("2 0 0 #nan #nan #nan\n", ost.str().c_str());
}

TEST(ErrorModel, error_model_read)
{
	//	 0 0 0 0 0 #nan
	//	1 0 0 0 #nan #nan
	//	2 0 0 #nan #nan #nan
	//	3 0 #nan #nan #nan #nan
	//	4 #nan #nan #nan #nan #nan
	std::istringstream ist("maxcnt:4\ncntdiff 1 2 3 4 5\n 0 0 0 0 0 #nan\n");
	ErrorStruct err;
	ist >> err;
	LONGS_EQUAL(4, err.maxfamilysize);
	LONGS_EQUAL(1, err.fromdiff);
	LONGS_EQUAL(5, err.todiff);
}

TEST(ErrorModel, add_and_remove_error_model)
{
	pCafeTree tree = create_tree();

	pCafeFamily pfamily = cafe_family_init({"chimp", "human", "mouse", "rat", "dog"});
	cafe_family_set_species_index(pfamily, tree);

	ErrorStruct err;
	pfamily->error_ptr = NULL;
	pErrorStruct perr = &err;
	init_error_ptr(pfamily, tree, &err, "human");
	CHECK(pfamily->error_ptr != NULL);
	POINTERS_EQUAL(perr, pfamily->error_ptr[1]);

	pfamily->errors = arraylist_new(10);
	arraylist_add(pfamily->errors, perr);
	remove_error_model(pfamily, tree, "human");
	POINTERS_EQUAL(NULL, pfamily->error_ptr[1]);

  cafe_family_free(pfamily);
}

