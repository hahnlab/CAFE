#include <sstream>
#include <vector>
#include <fstream>


#include "CppUTest/TestHarness.h"

#include <error_model.h>
#include <gene_family.h>
#include "cafe_commands.h"

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
  err.maxfamilysize = 0;
	ist >> err;
	LONGS_EQUAL(4, err.maxfamilysize);
	LONGS_EQUAL(1, err.fromdiff);
	LONGS_EQUAL(5, err.todiff);
}

TEST(ErrorModel, ErrorStruct_read_maxfamilysize_should_not_shrink)
{
  //	 0 0 0 0 0 #nan
  //	1 0 0 0 #nan #nan
  //	2 0 0 #nan #nan #nan
  //	3 0 #nan #nan #nan #nan
  //	4 #nan #nan #nan #nan #nan
  std::istringstream ist("maxcnt:4\ncntdiff 1 2 3 4 5\n 0 0 0 0 0 #nan\n");
  ErrorStruct err;
  err.maxfamilysize = 10;
  ist >> err;
  LONGS_EQUAL(10, err.maxfamilysize);
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

TEST(ErrorModel, read_freq_from_measures_errors)
{
	std::vector<int> sizeFreq(100);
	std::ifstream nil("nonexistent.tab");
	std::istringstream empty;
	const char *tab1 = "desc\tid\ta\tb\tc\nx\ty\t1\t2\t3\nz\tq\t4\t5\t6";
	const char *tab2 = "desc\tid\ta\tb\tc\nx\ty\t1\t2\t3\nz\tq\t4\t5\t3";
	std::istringstream t1(tab1);
	std::istringstream t2(tab2);
	CHECK_THROWS(io_error, read_freq_from_measures(&nil, &t2, &sizeFreq[0], 0));
	CHECK_THROWS(io_error, read_freq_from_measures(&t1, &nil, &sizeFreq[0], 0));

	std::istringstream t3(tab2);
	CHECK_THROWS(io_error, read_freq_from_measures(&empty, &t3, &sizeFreq[0], 0));
	std::istringstream t4(tab1);
	std::istringstream empty2;
	CHECK_THROWS(io_error, read_freq_from_measures(&t4, &empty2, &sizeFreq[0], 0));
}
	
TEST(ErrorModel, read_freq_from_measures)
{
	std::istringstream t1("desc\tid\ta\tb\tc\nx\ty\t1\t2\t3\nz\tq\t4\t5\t6");
	std::istringstream t2("desc\tid\ta\tb\tc\nx\ty\t1\t2\t3\nz\tq\t4\t5\t3");

	std::vector<int> sizeFreq(100);
	int max = read_freq_from_measures(&t1, &t2, &sizeFreq[0], 0);
	LONGS_EQUAL(2, sizeFreq[1]);
	LONGS_EQUAL(2, sizeFreq[2]);
	LONGS_EQUAL(3, sizeFreq[3]);
	LONGS_EQUAL(2, sizeFreq[4]);
	LONGS_EQUAL(2, sizeFreq[5]);
	LONGS_EQUAL(1, sizeFreq[6]);
	LONGS_EQUAL(6, max);
}

TEST(ErrorModel, read_freq_from_measures_arg_is_max)
{
	std::istringstream t1("desc\tid\ta\tb\tc\nx\ty\t1\t2\t3\nz\tq\t4\t5\t6");
	std::istringstream t2("desc\tid\ta\tb\tc\nx\ty\t1\t2\t3\nz\tq\t4\t5\t3");

	std::vector<int> sizeFreq(100);
	int max = read_freq_from_measures(&t1, &t2, &sizeFreq[0], 10);
	LONGS_EQUAL(10, max);
}

TEST(ErrorModel, read_freq_from_measures_second_file_optional)
{
	std::istringstream t2("desc\tid\ta\tb\tc\nx\ty\t1\t2\t3\nz\tq\t4\t5\t3");
	std::vector<int> sizeFreq(100);
	int max = read_freq_from_measures(&t2, NULL, &sizeFreq[0], 0);
	LONGS_EQUAL(1, sizeFreq[1]);
	LONGS_EQUAL(1, sizeFreq[2]);
	LONGS_EQUAL(2, sizeFreq[3]);
	LONGS_EQUAL(1, sizeFreq[4]);
	LONGS_EQUAL(1, sizeFreq[5]);
	LONGS_EQUAL(0, sizeFreq[6]);
	LONGS_EQUAL(5, max);
}

