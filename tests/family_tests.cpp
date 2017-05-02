#include <sstream>
#include <fstream>

#include "CppUTest/TestHarness.h"
#include "gene_family.h"

#include "cafe_commands.h"

extern "C" {
#include <utils_string.h>
#include <cafe_shell.h>
#include <tree.h>
#include <cafe.h>
};

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

TEST(FamilyTests, TestCafeFamily)
{
	std::ifstream ifst("Nonexistent.tab");
	POINTERS_EQUAL(NULL, load_gene_families(ifst, 1, '\t'));

  pCafeFamily fam = cafe_family_init({ "Dog", "Chimp", "Human", "Mouse", "Rat" });

  cafe_family_add_item(fam, { "OTOPETRIN", "ENSF1", "3", "5", "7", "11", "13" });

	LONGS_EQUAL(5, fam->num_species);
	LONGS_EQUAL(1, fam->flist->size);

	STRCMP_EQUAL("Dog", fam->species[0]);
	STRCMP_EQUAL("Rat", fam->species[4]);

	pCafeFamilyItem item = cafe_family_get_family_item(fam, "ENSF1");
	CHECK(item != NULL);
	LONGS_EQUAL(3, item->count[0]);
	LONGS_EQUAL(13, item->count[4]);
	STRCMP_EQUAL("OTOPETRIN", item->desc);

  cafe_family_free(fam);

}

TEST(FamilyTests, load_gene_families_tab_separated)
{
  std::ifstream ifst("Nonexistent.tab");
  POINTERS_EQUAL(NULL, load_gene_families(ifst, 1, '\t'));

  std::istringstream ist("FAMILYDESC\tFAMILY\tA\tB\nOTOPETRIN\tENSF2\t2\t10\n");
  pCafeFamily fam = load_gene_families(ist, 0, '\t');
  LONGS_EQUAL(2, fam->num_species);
  STRCMP_EQUAL("A", fam->species[0]);
  STRCMP_EQUAL("B", fam->species[1]);

  LONGS_EQUAL(1, fam->flist->size);

  pCafeFamilyItem item = cafe_family_get_family_item(fam, "ENSF2");
  CHECK(item != NULL);
  LONGS_EQUAL(2, item->count[0]);
  LONGS_EQUAL(10, item->count[1]);
  STRCMP_EQUAL("OTOPETRIN", item->desc);
  cafe_family_free(fam);
}

TEST(FamilyTests, load_gene_families_comma_separated)
{
  std::istringstream ist("FAMILY DESC,FAMILY,A,B\nOXYSTEROL BINDING RELATED,ENSF2,2,10\n");
  pCafeFamily fam = load_gene_families(ist, 0, ',');
  LONGS_EQUAL(2, fam->num_species);
  STRCMP_EQUAL("A", fam->species[0]);
  STRCMP_EQUAL("B", fam->species[1]);

  LONGS_EQUAL(1, fam->flist->size);

  pCafeFamilyItem item = cafe_family_get_family_item(fam, "ENSF2");
  CHECK(item != NULL);
  LONGS_EQUAL(2, item->count[0]);
  LONGS_EQUAL(10, item->count[1]);
  STRCMP_EQUAL("OXYSTEROL BINDING RELATED", item->desc);
  cafe_family_free(fam);
}

TEST(FamilyTests, load_gene_families_is_not_whitespace_separated)
{
  std::istringstream ist("FAMILYDESC FAMILY A B C\n");
  try
  {
    pCafeFamily fam = load_gene_families(ist, 0, '\t');
    FAIL("Exception should have been thrown");
    cafe_family_free(fam);
  }
  catch (std::runtime_error& err)
  {
    STRCMP_EQUAL("Failed to identify species for gene families", err.what());
  }
}

TEST(FamilyTests, cafe_family_set_size)
{
	pCafeFamily pcf = cafe_family_init({"chimp", "human", "mouse", "rat", "dog" });
	cafe_family_add_item(pcf, { "description", "id", "3", "5", "7", "11", "13" });

	pCafeTree ptree = create_tree();
	cafe_family_set_species_index(pcf, ptree);
	cafe_family_set_size(pcf, 0, ptree);

	// leaf nodes are at 0,2,4,6,8. They should match the sizes given above
	LONGS_EQUAL(3, ((pCafeNode)ptree->super.nlist->array[0])->familysize);
	LONGS_EQUAL(5, ((pCafeNode)ptree->super.nlist->array[2])->familysize);
	LONGS_EQUAL(7, ((pCafeNode)ptree->super.nlist->array[4])->familysize);
	LONGS_EQUAL(11, ((pCafeNode)ptree->super.nlist->array[6])->familysize);
	LONGS_EQUAL(13, ((pCafeNode)ptree->super.nlist->array[8])->familysize);

  cafe_family_free(pcf);
};

TEST(FamilyTests, cafe_family_init)
{
	pCafeFamily pcf = cafe_family_init({"chimp", "human", "mouse", "rat", "dog" });
	LONGS_EQUAL(5, pcf->num_species);
	STRCMP_EQUAL("dog", pcf->species[4]);
	STRCMP_EQUAL("rat", pcf->species[3]);
	LONGS_EQUAL(5, pcf->num_species);
	LONGS_EQUAL(0, pcf->max_size);
	LONGS_EQUAL(0, pcf->flist->size);
	LONGS_EQUAL(-1, pcf->index[1]);
	LONGS_EQUAL(-1, pcf->index[2]);
  cafe_family_free(pcf);

}

TEST(FamilyTests, cafe_family_add_item)
{
	pCafeFamily pcf = cafe_family_init({"chimp", "human", "mouse", "rat", "dog" });

	cafe_family_add_item(pcf, { "description", "id", "3", "5", "7", "11", "13" });
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
  cafe_family_free(pcf);

}

TEST(FamilyTests, cafe_family_set_species_index)
{
	pCafeTree tree = create_tree();
	pCafeFamily pcf = cafe_family_init({"chimp", "human", "mouse", "rat", "dog" });

	LONGS_EQUAL(-1, pcf->index[0]);
	int errcode = cafe_family_set_species_index(pcf, tree);
	LONGS_EQUAL(0, errcode);
	LONGS_EQUAL(0, pcf->index[0]);
	LONGS_EQUAL(2, pcf->index[1]);
	LONGS_EQUAL(4, pcf->index[2]);
	LONGS_EQUAL(6, pcf->index[3]);
	LONGS_EQUAL(8, pcf->index[4]);

  cafe_family_free(pcf);

  pcf = cafe_family_init({ "chimp", "human", "mouse", "rat" });
	errcode = cafe_family_set_species_index(pcf, tree);
	LONGS_EQUAL(-1, errcode);	// Error because dog is missing from family list
  cafe_family_free(pcf);
}

TEST(FamilyTests, cafe_family_set_size_with_family_forced)
{
	// setup
	pCafeFamily pcf = cafe_family_init({"chimp", "human", "mouse", "rat", "dog" });
	pCafeTree cafe_tree = create_tree();
	cafe_family_set_species_index(pcf, cafe_tree);
	cafe_family_add_item(pcf, { "description", "id", "3", "5", "7", "11", "13" });

	// SUT
	cafe_family_set_size_with_family_forced(pcf, 0, cafe_tree);
	LONGS_EQUAL(1, cafe_tree->rootfamilysizes[0]);
	LONGS_EQUAL(16, cafe_tree->rootfamilysizes[1]);
	LONGS_EQUAL(63, cafe_tree->familysizes[1]);
	LONGS_EQUAL(16, cafe_tree->rfsize);
  cafe_family_free(pcf);
}

TEST(FamilyTests, write_family_gainloss)
{
	std::ostringstream ost;
	pCafeFamily pcf = cafe_family_init({"chimp", "human", "mouse", "rat", "dog" });
	cafe_family_add_item(pcf, { "description", "id", "3", "5", "7", "11", "13" });
	pCafeTree tree1 = create_tree();
	cafe_family_set_species_index(pcf, tree1);
	pCafeTree tree2 = cafe_tree_copy(tree1);
	int **viterbiNodeFamilysizes = (int**)memory_new_2dim(tree1->super.nlist->size, pcf->num_species, sizeof(int));
	cafe_family_set_size(pcf, 0, tree1);
	set_node_familysize(tree1, viterbiNodeFamilysizes, 0);
	int actual = write_family_gainloss(ost, "family_id", tree1, tree2);
	LONGS_EQUAL(39, actual);
	STRCMP_EQUAL("family_id\t39\t(((chimp_3<3>:6,human_5<5>:6)_0<0>:81,(mouse_7<7>:17,rat_11<11>:17)_0<0>:70)_0<0>:6,dog_13<13>:9)_0\n", ost.str().c_str())
    cafe_family_free(pcf);
}

TEST(FamilyTests, write_family)
{
	std::ostringstream ost;
	pCafeFamily pcf = cafe_family_init({"chimp", "human", "mouse", "rat", "dog" });
	cafe_family_add_item(pcf, { "my_description", "my_id", "3", "5", "7", "11", "13" });

	write_family(ost, pcf);
	STRCMP_CONTAINS("Desc\tFamily ID\tchimp\thuman\tmouse\trat\tdog\n", ost.str().c_str())
	STRCMP_CONTAINS("my_description\tmy_id\t3\t5\t7\t11\t13\n", ost.str().c_str())
  cafe_family_free(pcf);

}

