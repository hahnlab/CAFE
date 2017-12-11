#ifndef REPORTS_C6B4B223_5C6D_459E_B2D8_CAC46A40A5B4
#define REPORTS_C6B4B223_5C6D_459E_B2D8_CAC46A40A5B4

#include <vector>
#include <string>

extern "C" {
#include "cafe_shell.h"
#include "family.h"
}

#include "viterbi.h"


class Globals;
class viterbi_parameters;

struct family_line_item
{
	std::string node_id;
	std::string tree;
	double max_p_value;
	std::vector<std::pair<double, double> > pvalues;
	std::vector<double> cut_pvalues;
	std::vector<double> likelihood_ratios;

	family_line_item(pCafeFamily family, pCafeTree pcafe, double** likelihoodRatios, viterbi_parameters& viterbi, int i, std::string node_id);
	family_line_item() : max_p_value(0.0) {}
};

class tree_visualization
{
    virtual void serialize(std::ostream& ost) const = 0;
    friend std::ostream& operator<<(std::ostream& ost, const tree_visualization& newick);
};

class newick_visualization : public tree_visualization
{
    pTree _tree;
public:
    newick_visualization(pTree tree) : _tree(tree)
    {

    }
    virtual void serialize(std::ostream& ost) const;

};

struct coord
{
    double x;
    double y;
    coord() {}
    coord(double _x, double _y) : x(_x), y(_y)
    {

    }
};

class svg_visualization : public tree_visualization
{
    pTree _tree;
    int width;
    int left_margin;
    int right_margin;
    int top_margin;
    double legend_ratio;
    double _font_size;
    double precision;
    int tip_space;
    int longest_label;
    double scaler;
    void set_xcoord(pTreeNode node);
    void set_ycoord(pTreeNode node);
public:
    svg_visualization(pTree tree);
    virtual void serialize(std::ostream& ost) const;
    std::map<int, coord> coordinates;
    void plot_node(std::ostream& ost, pTreeNode node) const;
    int height;
    bool legend;

};

class ascii_visualization : public tree_visualization
{
    pTree _tree;
    int _width;
public:
    ascii_visualization(pTree tree, int width) : _tree(tree), _width(width) {}
    virtual void serialize(std::ostream& ost) const;
};

struct Report
{
  enum Formats { Unknown, Text, HTML, JSON };
  static int report_format;

  pCafeTree aTree;
	std::string tree;
	std::string lambda_tree;
	std::vector<double> lambdas;
	std::vector<double> averageExpansion;
	std::vector<change> changes;
	std::vector<std::pair<int, int> > node_pairs;
	std::vector<family_line_item> family_line_items;
	std::vector<int> branch_cutting_output_format;

	Report() {}
	Report(pCafeParam param, viterbi_parameters& viterbi);
};

struct report_parameters
{
  bool branchcutting;
  bool likelihood;
  bool lh2;
  bool just_save;
  Report::Formats format;
  std::string name;
};

void update_depths(pTreeNode node, std::map<int, double>& depths, double curr_depth);

report_parameters get_report_parameters(std::vector<std::string> tokens);
int cafe_cmd_report(Globals& globals, std::vector<std::string> tokens);
void write_viterbi(std::ostream& ost, const Report& viterbi);
void write_families_header(std::ostream& ost, bool cutPvalues, bool likelihoodRatios);
void cafe_do_report(Globals& globals, viterbi_parameters& viterbi, report_parameters* params);
int cafe_report_retrieve_data(const char* file, pCafeParam param, viterbi_parameters& viterbi);

std::ostream& operator<<(std::ostream& ost, const family_line_item& item);
std::ostream& operator<<(std::ostream& ost, const Report& report);

/// These I/O manipulators allow selecting the report format
std::ios_base& json(std::ios_base& os);
std::ios_base& html(std::ios_base& os);
#endif
