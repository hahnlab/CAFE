#include <vector>
#include <sstream>
#include <iostream>

#include "CppUTest/TestHarness.h"
#include "CppUTest/CommandLineTestRunner.h"

#include <reports.h>
#include <Globals.h>
#include <gene_family.h>

extern "C" {
//#include <utils_string.h>
//#include <cafe_shell.h>
//#include <tree.h>
#include <cafe.h>
//#include <chooseln_cache.h>

  extern pBirthDeathCacheArray probability_cache;
};

static pCafeTree create_tree(family_size_range range)
{
  const char *newick_tree = "(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)";
  char tree[100];
  strcpy(tree, newick_tree);
  pCafeTree pcafe = cafe_tree_new(tree, &range, 0.01, 0);

  return pcafe;
}

TEST_GROUP(ReportTests)
{
  void setup()
  {
    srand(10);
  }
};

TEST(ReportTests, get_report_parameters)
{
  std::vector<std::string> tokens;
  tokens.push_back("report");
  tokens.push_back("myreport");
  tokens.push_back("branchcutting");

  report_parameters params = get_report_parameters(tokens);
  STRCMP_EQUAL("myreport", params.name.c_str());
  CHECK(params.branchcutting);
  CHECK_FALSE(params.likelihood);
  LONGS_EQUAL(Report::Text, params.format);

  tokens[2] = "lh2";
  params = get_report_parameters(tokens);
  CHECK(params.lh2);
  CHECK_FALSE(params.likelihood);

  tokens.push_back("html");
  params = get_report_parameters(tokens);
  LONGS_EQUAL(Report::HTML, params.format);

  tokens[3] = "json";
  params = get_report_parameters(tokens);
  LONGS_EQUAL(Report::JSON, params.format);
}

TEST(ReportTests, write_report)
{
  Globals globals;
  std::ostringstream ost;
  CafeFamily fam;
  fam.flist = arraylist_new(1);

  globals.viterbi->num_nodes = 0;
  family_size_range range;
  range.min = 0;
  range.max = 15;
  range.root_min = 0;
  range.root_max = 15;

  globals.param.pcafe = create_tree(range);
  globals.param.pfamily = &fam;
  probability_cache = NULL;
  globals.param.num_lambdas = 3;
  double lambdas[] = { 1.5, 2.5, 3.5 };
  globals.param.lambda = lambdas;
  globals.param.lambda_tree = (pTree)globals.param.pcafe;
  Report r(&globals.param, *globals.viterbi);
  ost << r;
  STRCMP_CONTAINS("Tree:(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)\n", ost.str().c_str());
  STRCMP_CONTAINS("Lambda:\t1.5\t2.5\t3.5\n", ost.str().c_str());
  STRCMP_CONTAINS("Lambda tree:\t(((:6,:6):81,(:17,:17):70):6,:9)\n", ost.str().c_str());
  STRCMP_CONTAINS("IDs of nodes:(((chimp<0>,human<2>)<1>,(mouse<4>,rat<6>)<5>)<3>,dog<8>)<7>\n", ost.str().c_str());
  STRCMP_CONTAINS("# Output format for: ' Average Expansion', 'Expansions', 'No Change', 'Contractions', and 'Branch-specific P-values' = (node ID, node ID): ", ost.str().c_str());
  STRCMP_CONTAINS("(0,2) (1,5) (4,6) (3,8) \n", ost.str().c_str());
  STRCMP_CONTAINS("# Output format for 'Branch cutting P-values' and 'Likelihood Ratio Test': (0, 1, 2, 3, 4, 5, 6, 7, 8)\n", ost.str().c_str());
  STRCMP_CONTAINS("", ost.str().c_str());
}

Report create_report()
{
  Report r;
  r.tree = "(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)";
  r.lambdas = { 0.1, 0.3, 0.5, 0.7 };
  r.lambda_tree = "(((:6,:6):81,(:17,:17):70):6,:9)";
  r.id_tree = "(((chimp<0>,human<2>)<1>,(mouse<4>,rat<6>)<5>)<3>,dog<8>)<7>";
  r.node_pairs = { { 1,2 },{ 3,4 },{ 5,6 } };
  r.averageExpansion = { 1.1, 1.3, 1.5, 1.7, 2.1, 2.3, 2.5, 2.7 };
  r.changes = { change(1,2,3), change(4,5,6), change(7,8,9), change(10,11,12) };
  family_line_item f;
  f.node_id = "C1";
  f.tree = "tree";
  f.max_p_value = 0.08;
  f.pvalues = { { 11,12 },{ 13,14 },{ 15,16 } };
  r.family_line_items = { f };

  return r;
}

TEST(ReportTests, write_html_report)
{
  Report r = create_report();

  std::ostringstream ost;

  ost << html << r;
  STRCMP_CONTAINS("<h3>Input tree</h3><p>(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)</p>", ost.str().c_str());
  STRCMP_CONTAINS("<h3>Lambda values</h3><ol><li>0.1</li><li>0.3</li><li>0.5</li><li>0.7</li></ol>", ost.str().c_str());
  STRCMP_CONTAINS("<h3>Lambda tree</h3><p>(((:6,:6):81,(:17,:17):70):6,:9)</p>", ost.str().c_str());
  STRCMP_CONTAINS("<tr><th>Node Pairs</th><th>(1,2)</th><th>(3,4)</th><th>(5,6)</th></tr>", ost.str().c_str());
  STRCMP_CONTAINS("<tr><td>Average Expansion</td><td>(1.1,1.3)</td><td>(1.5,1.7)</td><td>(2.1,2.3)</td><td>(2.5,2.7)</td></tr>", ost.str().c_str());
  STRCMP_CONTAINS("<tr><td>Expanded</td><td>(1,4)</td><td>(7,10)</td></tr>", ost.str().c_str());
  STRCMP_CONTAINS("<tr><td>Decreased</td><td>(3,6)</td><td>(9,12)</td></tr>", ost.str().c_str());
  STRCMP_CONTAINS("<h3>Families</h3><table><tr><th>ID</th><th>Newick</th><th>Family-wide P-value</th><th>Viterbi P-values</th><th>cut P-value</th><th>Likelihood Ratio</th></tr>", ost.str().c_str());
  STRCMP_CONTAINS("<tr><td>C1</td><td>tree</td><td>0.08</td><td>(11,12),(13,14),(15,16)</td></tr>", ost.str().c_str());
  //STRCMP_CONTAINS("\"Expanded\":[[1.1,1.3],[1.5,1.7],[2.1,2.3],[2.5,2.7]]", ost.str().c_str());

}

TEST(ReportTests, write_json_report)
{
  Report r = create_report();

  std::ostringstream ost;

  ost << json << r;
  STRCMP_CONTAINS("\"InputTree\":\"(((chimp:6,human:6):81,(mouse:17,rat:17):70):6,dog:9)\"", ost.str().c_str());
  STRCMP_CONTAINS("\"Lambdas\":[0.1,0.3,0.5,0.7]", ost.str().c_str());
  STRCMP_CONTAINS("\"LambdaTree\":\"(((:6,:6):81,(:17,:17):70):6,:9)\"", ost.str().c_str());
  STRCMP_CONTAINS("\"AverageExpansion\":[[1.1,1.3],[1.5,1.7],[2.1,2.3],[2.5,2.7]]", ost.str().c_str());
  STRCMP_CONTAINS("\"Expanded\":[[1,4],[7,10]]", ost.str().c_str());
  STRCMP_CONTAINS("\"Unchanged\":[[2,5],[8,11]]", ost.str().c_str());
  STRCMP_CONTAINS("\"Decreased\":[[3,6],[9,12]]", ost.str().c_str());
  STRCMP_CONTAINS("\"Families\":[{ \"Node\":\"C1\",\n\"Tree\":\"tree\",\n\"MaxPValue\":0.08,\n", ost.str().c_str());
  STRCMP_CONTAINS("\"PValues\":[[11,12],[13,14],[15,16]]", ost.str().c_str());
  //STRCMP_CONTAINS("\"Expanded\":[[1.1,1.3],[1.5,1.7],[2.1,2.3],[2.5,2.7]]", ost.str().c_str());

}

TEST(ReportTests, write_viterbi)
{
  std::ostringstream ost;
  Report r;
  for (int i = 1; i < 7; ++i)
    r.averageExpansion.push_back(i + 0.3);

  r.changes.push_back(change(1, 17, 41));
  r.changes.push_back(change(3, 19, 43));
  r.changes.push_back(change(5, 23, 47));
  r.changes.push_back(change(7, 29, 53));
  r.changes.push_back(change(11, 31, 59));
  r.changes.push_back(change(13, 37, 61));

  write_viterbi(ost, r);

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
  Globals globals;
  globals.param.likelihoodRatios = NULL;
  family_size_range range;
  range.min = range.root_min = 0;
  range.max = range.root_max = 10;
  pCafeTree tree = create_tree(range);
  globals.param.pcafe = tree;
  globals.param.pfamily = cafe_family_init({ "chimp", "human", "mouse", "rat", "dog" });
  cafe_family_set_species_index(globals.param.pfamily, tree);
  cafe_family_add_item(globals.param.pfamily, { "description", "id", "3", "5", "7", "11", "13" });

  viterbi_parameters* v = globals.viterbi;
  v->num_nodes = 6;

  pCafeNode pnode = (pCafeNode)tree->super.nlist->array[0];
  pCafeFamilyItem pitem = (pCafeFamilyItem)globals.param.pfamily->flist->array[0];
  v->viterbiPvalues[viterbi_parameters::NodeFamilyKey(pnode, pitem)] = .025;
  v->cutPvalues = NULL;

  double maxP = .1;
  v->maximumPvalues = &maxP;

  std::ostringstream ost;
  family_line_item item(globals.param.pfamily, globals.param.pcafe, globals.param.likelihoodRatios, *globals.viterbi, 0, "NodeZero");
  ost << item;
  STRCMP_EQUAL("NodeZero\t(((chimp_3:6,human_5:6)_0:81,(mouse_7:17,rat_11:17)_0:70)_0:6,dog_13:9)_0\t0.1\t((0.025,0),(0,0),(0,0))\t", ost.str().c_str());

  v->cutPvalues = (double**)memory_new_2dim(6, 1, sizeof(double));
  std::ostringstream ost2;
  family_line_item item2(globals.param.pfamily, globals.param.pcafe, globals.param.likelihoodRatios, *globals.viterbi, 0, "NodeZero");
  ost2 << item2;
  STRCMP_EQUAL("NodeZero\t(((chimp_3:6,human_5:6)_0:81,(mouse_7:17,rat_11:17)_0:70)_0:6,dog_13:9)_0\t0.1\t((0.025,0),(0,0),(0,0))\t(0,0,0,0,0,0)\t", ost2.str().c_str());

  v->cutPvalues = NULL;
  globals.param.likelihoodRatios = (double**)memory_new_2dim(tree->super.size, 1, sizeof(double));
  std::ostringstream ost3;
  family_line_item item3(globals.param.pfamily, globals.param.pcafe, globals.param.likelihoodRatios, *globals.viterbi, 0, "NodeZero");
  ost3 << item3;
  STRCMP_EQUAL("NodeZero\t(((chimp_3:6,human_5:6)_0:81,(mouse_7:17,rat_11:17)_0:70)_0:6,dog_13:9)_0\t0.1\t((0.025,0),(0,0),(0,0))\t(0,0,0,0,0,0,0,0,0)", ost3.str().c_str());

  cafe_family_free(globals.param.pfamily);
}

