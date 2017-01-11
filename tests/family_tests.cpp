#include <sstream>
#include "CppUTest/TestHarness.h"

#include "cafe_commands.h"

extern "C" {
#include <utils_string.h>
#include <cafe_shell.h>
#include <tree.h>
#include <cafe.h>
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


static pCafeTree create_tree()
{
	const char *newick_tree = "(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)";
	char tree[100];
	strcpy(tree, newick_tree);
	family_size_range range;
	range.min = range.root_min = 1;
	range.max = range.root_max = 1;

	return cafe_tree_new(tree, &range, 0, 0);
}


TEST_GROUP(FamilyTests)
{
};

TEST(FamilyTests, TestCafeFamilyNew)
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


TEST(FamilyTests, cafe_family_set_size)
{
	const char *species[] = { "", "", "chimp", "human", "mouse", "rat", "dog" };
	pCafeFamily pcf = cafe_family_init(build_arraylist(species, 7));
	const char *values[] = { "description", "id", "3", "5", "7", "11", "13" };
	cafe_family_add_item(pcf, build_arraylist(values, 7));

	pCafeTree ptree = create_tree();
	cafe_family_set_species_index(pcf, ptree);
	cafe_family_set_size(pcf, 0, ptree);

	// leaf nodes are at 0,2,4,6,8. They should match the sizes given above
	LONGS_EQUAL(3, ((pCafeNode)ptree->super.nlist->array[0])->familysize);
	LONGS_EQUAL(5, ((pCafeNode)ptree->super.nlist->array[2])->familysize);
	LONGS_EQUAL(7, ((pCafeNode)ptree->super.nlist->array[4])->familysize);
	LONGS_EQUAL(11, ((pCafeNode)ptree->super.nlist->array[6])->familysize);
	LONGS_EQUAL(13, ((pCafeNode)ptree->super.nlist->array[8])->familysize);
};

TEST(FamilyTests, cafe_family_init)
{
	const char *species[] = { "", "", "chimp", "human", "mouse", "rat", "dog" };
	pArrayList psplit = build_arraylist(species, 7);

	pCafeFamily pcf = cafe_family_init(psplit);
	arraylist_free(psplit, NULL);
	LONGS_EQUAL(5, pcf->num_species);
	STRCMP_EQUAL("dog", pcf->species[4]);
	STRCMP_EQUAL("rat", pcf->species[3]);
	LONGS_EQUAL(5, pcf->num_species);
	LONGS_EQUAL(0, pcf->max_size);
	LONGS_EQUAL(0, pcf->flist->size);
	LONGS_EQUAL(-1, pcf->index[1]);
	LONGS_EQUAL(-1, pcf->index[2]);

}

TEST(FamilyTests, cafe_family_add_item)
{
	const char *species[] = { "", "", "chimp", "human", "mouse", "rat", "dog" };

	pArrayList psplit = build_arraylist(species, 7);
	pCafeFamily pcf = cafe_family_init(psplit);

	const char *values[] = { "description", "id", "3", "5", "7", "11", "13" };
	cafe_family_add_item(pcf, build_arraylist(values, 7));
	arraylist_free(psplit, NULL);
	LONGS_EQUAL(1, pcf->flist->size);
	pCafeFamilyItem pitem = (pCafeFamilyItem)arraylist_get(pcf->flist, 0);
	STRCMP_EQUAL("description", pitem->desc);
	STRCMP_EQUAL("id", pitem->id);
	LONGS_EQUAL(3, pitem->count[0]);
	LONGS_EQUAL(5, pitem->count[1]);
	LONGS_EQUAL(7, pitem->count[2]);
	LONGS_EQUAL(11, pitem->count[3]);
	LONGS_EQUAL(13, pitem->count[4]);
	LONGS_EQUAL(pitem->maxlh, -1);
	LONGS_EQUAL(pitem->ref, -1);
	LONGS_EQUAL(pitem->lambda, NULL);
	LONGS_EQUAL(pitem->mu, NULL);
	LONGS_EQUAL(pitem->holder, 1);
	LONGS_EQUAL(13, pcf->max_size);

}

TEST(FamilyTests, cafe_family_set_species_index)
{
	pCafeTree tree = create_tree();
	const char *species[] = { "", "", "chimp", "human", "mouse", "rat", "dog" };
	pCafeFamily pcf = cafe_family_init(build_arraylist(species, 7));

	LONGS_EQUAL(-1, pcf->index[0]);
	int errcode = cafe_family_set_species_index(pcf, tree);
	LONGS_EQUAL(0, errcode);
	LONGS_EQUAL(0, pcf->index[0]);
	LONGS_EQUAL(2, pcf->index[1]);
	LONGS_EQUAL(4, pcf->index[2]);
	LONGS_EQUAL(6, pcf->index[3]);
	LONGS_EQUAL(8, pcf->index[4]);

	pcf = cafe_family_init(build_arraylist(species, 6));
	errcode = cafe_family_set_species_index(pcf, tree);
	LONGS_EQUAL(-1, errcode);	// Error because dog is missing from family list
}

TEST(FamilyTests, cafe_family_set_size_with_family_forced)
{
	// setup
	const char *species[] = { "", "", "chimp", "human", "mouse", "rat", "dog" };
	pCafeFamily pcf = cafe_family_init(build_arraylist(species, 7));
	pCafeTree cafe_tree = create_tree();
	cafe_family_set_species_index(pcf, cafe_tree);
	const char *values[] = { "description", "id", "3", "5", "7", "11", "13" };
	cafe_family_add_item(pcf, build_arraylist(values, 7));

	// SUT
	cafe_family_set_size_with_family_forced(pcf, 0, cafe_tree);
	LONGS_EQUAL(1, cafe_tree->rootfamilysizes[0]);
	LONGS_EQUAL(16, cafe_tree->rootfamilysizes[1]);
	LONGS_EQUAL(63, cafe_tree->familysizes[1]);
	LONGS_EQUAL(16, cafe_tree->rfsize);
}

TEST(FamilyTests, write_family_gainloss)
{
	std::ostringstream ost;
	const char *species[] = { "", "", "chimp", "human", "mouse", "rat", "dog" };
	pCafeFamily pcf = cafe_family_init(build_arraylist(species, 7));
	const char *values[] = { "description", "id", "3", "5", "7", "11", "13" };
	cafe_family_add_item(pcf, build_arraylist(values, 7));
	pCafeTree tree1 = create_tree();
	cafe_family_set_species_index(pcf, tree1);
	pCafeTree tree2 = cafe_tree_copy(tree1);
	int **viterbiNodeFamilysizes = (int**)memory_new_2dim(tree1->super.nlist->size, pcf->num_species, sizeof(int));
	cafe_family_set_size(pcf, 0, tree1);
	set_node_familysize(tree1, viterbiNodeFamilysizes, 0);
	int actual = write_family_gainloss(ost, "family_id", tree1, tree2);
	LONGS_EQUAL(39, actual);
	STRCMP_EQUAL("family_id\t39\t(((chimp_3<3>:6,human_5<5>:6)_0<0>:81,(mouse_7<7>:17,rat_11<11>:17)_0<0>:70)_0<0>:6,dog_13<13>:9)_0\n", ost.str().c_str())
}

TEST(FamilyTests, write_family)
{
	std::ostringstream ost;
	const char *species[] = { "", "", "chimp", "human", "mouse", "rat", "dog" };
	pCafeFamily pcf = cafe_family_init(build_arraylist(species, 7));
	const char *values[] = { "my_description", "my_id", "3", "5", "7", "11", "13" };
	cafe_family_add_item(pcf, build_arraylist(values, 7));

	write_family(ost, pcf);
	STRCMP_CONTAINS("Desc\tFamily ID\tchimp\thuman\tmouse\trat\tdog\n", ost.str().c_str())
	STRCMP_CONTAINS("my_description\tmy_id\t3\t5\t7\t11\t13\n", ost.str().c_str())
}

