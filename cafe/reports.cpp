#include <stdexcept>
#include <ostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <algorithm>

#include <strings.h>

#include "reports.h"
#include "likelihood_ratio.h"
#include "pvalue.h"
#include "branch_cutting.h"
#include "viterbi.h"
#include "Globals.h"

using namespace std;

extern "C" {
#include "cafe.h"

	extern pTree tmp_lambda_tree;
	void cafe_shell_set_lambda(pCafeParam param, double* lambda);
	size_t file_read_line(pString pstr, FILE* fp);
	void cafe_report_load_bc_or_lhr_list(char* data, double** pvalues, int i, int nnodes);
	void phylogeny_lambda_parse_func(pTree ptree, pTreeNode ptnode);
	double* cafe_report_load_data_double_list(char* data, int delimiter, int* num);
	double* cafe_report_load_data_double_pairs(char* data, int delimiter);
	int* cafe_report_load_data_int_pairs(char* data, int delimiter);
}

/// Writes a pair of arbitrary objects to a stream, surrounded by brackets in JSON mode 
/// or parenthesis otherwise, and separated by a comma
template<typename T, typename U>
std::ostream& operator<<(ostream& ost, const std::pair<T, U>& pair)
{
  Report::Formats format = static_cast<Report::Formats>(ost.iword(Report::report_format));
  if (format == Report::JSON)
    ost << "[" << pair.first << "," << pair.second << "]";
  else
    ost << "(" << pair.first << "," << pair.second << ")";
  return ost;
}

/// Writes a vector of arbitrary objects to a stream, separated by commas
template<typename T>
std::ostream& operator<<(ostream& ost, const std::vector<T>& v)
{
  for (size_t i = 0; i < v.size() - 1; i++)
  {
    ost << v[i] << ",";
  }
  ost << v[v.size() - 1];

  return ost;
}

/// Converts a vector of objects to pairs of objects
template<typename T>
std::vector<pair<T, T> > to_pairs(const std::vector<T>& v)
{
  vector<pair<T, T> > x;
  for (size_t b = 0; b < v.size() / 2; b++)
  {
    x.push_back(pair<T, T>(v[2 * b], v[2 * b + 1]));
  }
  return x;
}

/// Accepts a change structure and returns its expand value
int expanded(change c)
{
  return c.expand;
}

/// Accepts a change structure and returns its unchanged value
int unchanged(change c)
{
  return c.remain;
}

/// Accepts a change structure and returns its decrease value
int decreased(change c)
{
  return c.decrease;
}

template<typename func>
std::vector<int> get_change(const std::vector<change>& v, func f)
{
  vector<int> x(v.size());
  transform(v.begin(), v.end(), x.begin(), f);
  return x;
}


void lambda_tree_string(pString pstr, pPhylogenyNode pnode)
{
	if (pnode->taxaid != -1)
	{
		string_fadd(pstr, "%d", pnode->taxaid + 1);
	}
}

void cafe_report_set_viterbi(pCafeFamily family, pCafeTree pcafe, viterbi_parameters& viterbi, int i)
{
	cafe_family_set_size(family, i, pcafe);
  pCafeFamilyItem pitem = (pCafeFamilyItem)family->flist->array[i];
  viterbi.set_node_familysize(pcafe, pitem);
}

void cafe_tree_string_id(pString pstr, pPhylogenyNode pnode)
{
	if (pnode->name) string_fadd(pstr, "%s", pnode->name);
	string_fadd(pstr, "<%d>", pnode->super.id);
}

string quoted(string s)
{
  return "\"" + s + "\"";
}

void write_viterbi(ostream& ost, const Report& viterbi)
{
  vector<pair<double, double> > avg_expansion = to_pairs(viterbi.averageExpansion);
  vector<pair<int, int> > expand = to_pairs(get_change(viterbi.changes, expanded));
  vector<pair<int, int> > remain = to_pairs(get_change(viterbi.changes, unchanged));
  vector<pair<int, int> > decrease = to_pairs(get_change(viterbi.changes, decreased));
  Report::Formats format = static_cast<Report::Formats>(ost.iword(Report::report_format));
  switch (format)
  {
  case Report::Text:
  case Report::Unknown:
    ost << "Average Expansion:";
    for (size_t b = 0; b < avg_expansion.size(); b++)
    {
      ost << "\t" << avg_expansion[b];
    }
    ost << "\n";
    ost << "Expansion :";
    for (size_t b = 0; b < expand.size(); b++)
    {
      ost << "\t" << expand[b];
    }
    ost << "\n";

    ost << "nRemain :";
    for (size_t b = 0; b < remain.size(); b++)
    {
      ost << "\t" << remain[b];
    }
    ost << "\n";

    ost << "nDecrease :";
    for (size_t b = 0; b < decrease.size(); b++)
    {
      ost << "\t" << decrease[b];
    }
    ost << "\n";
    break;
  case Report::HTML:
    ost << "<td>Average Expansion</td>";
    for (size_t b = 0; b < avg_expansion.size(); b++)
    {
      ost << "<td>" << avg_expansion[b] << "</td>";
    }
    ost << "</tr><tr>";
    ost << "<td>Expanded</td>";
    for (size_t b = 0; b < expand.size(); b++)
    {
      ost << "<td>" << expand[b] << "</td>";
    }
    ost << "</tr><tr>";
    ost << "<td>Remained Same</td>";
    for (size_t b = 0; b < remain.size(); b++)
    {
      ost << "<td>" << remain[b] << "</td>";
    }
    ost << "</tr><tr>";
    ost << "<td>Decreased</td>";
    for (size_t b = 0; b < decrease.size(); b++)
    {
      ost << "<td>" << decrease[b] << "</td>";
    }
    break;
  case Report::JSON:
    ost << quoted("AverageExpansion") << ":[" << avg_expansion << "]," << endl;
    ost << quoted("Expanded") << ":[" << to_pairs(get_change(viterbi.changes, expanded)) << "]," << endl;
    ost << quoted("Unchanged") << ":[" << to_pairs(get_change(viterbi.changes, unchanged)) << "]," << endl;
    ost << quoted("Decreased") << ":[" << to_pairs(get_change(viterbi.changes, decreased)) << "]," << endl;
    break;

  }
}

void write_families_header(ostream& ost, bool cutPvalues, bool likelihoodRatios)
{
	ost << "'ID'\t'Newick'\t'Family-wide P-value'\t'Viterbi P-values'";
	if (cutPvalues)
		ost << "\t'cut P-value'";
	if (likelihoodRatios)
		ost << "\t'Likelihood Ratio'";
	ost << "\n";
}

void write_doubles(ostream&ost, vector<double> items)
{
	ost << "(";

	vector<double>::iterator iter = items.begin();
	while (iter != items.end())
	{
		if (*iter == -1)
		{
			ost << "-";
		}
		else
		{
			ost << *iter;
		}
		*iter++;
		if (iter != items.end())
			ost << ",";
	}
	ost << ")";
}

family_line_item::family_line_item(pCafeFamily family, pCafeTree pcafe, double** likelihoodRatios, viterbi_parameters& viterbi, int i, string node_id)
{
	this->node_id = node_id;
	cafe_report_set_viterbi(family, pcafe, viterbi, i);
	pString pstr = cafe_tree_string(pcafe);
	tree = pstr->buf;
	string_free(pstr);
	max_p_value = viterbi.maximumPvalues[i];
	for (int b = 0; b < viterbi.num_nodes / 2; b++)
	{
    pCafeNode node1 = (pCafeNode)pcafe->super.nlist->array[2 * b];
    pCafeNode node2 = (pCafeNode)pcafe->super.nlist->array[2 * b +1];
    pCafeFamilyItem item = (pCafeFamilyItem)family->flist->array[i];
    double v1 = viterbi.viterbiPvalues[viterbi_parameters::NodeFamilyKey(node1, item)];
    double v2 = viterbi.viterbiPvalues[viterbi_parameters::NodeFamilyKey(node2, item)];
    pvalues.push_back(std::pair<double, double>(v1, v2));
	}

	if (viterbi.cutPvalues)
	{
		cut_pvalues.resize(viterbi.num_nodes);
		for (int b = 0; b < viterbi.num_nodes; b++)
		{
			cut_pvalues[b] = viterbi.cutPvalues[b][i];
		}
	}

	if (likelihoodRatios)
	{
		likelihood_ratios.resize(pcafe->super.nlist->size);
		for (size_t b = 0; b < likelihood_ratios.size(); b++)
		{
			likelihood_ratios[b] = likelihoodRatios[b][i];
		}
	}
}

// allocates the iword storage for use with Report objects
int Report::report_format = std::ios_base::xalloc();

// This I/O manipulator to select format types
std::ios_base& json(std::ios_base& os)
{
  os.iword(Report::report_format) = Report::JSON;
  return os;
}

std::ios_base& html(std::ios_base& os)
{
  os.iword(Report::report_format) = Report::HTML;
  return os;
}

std::ostream& operator<<(std::ostream& ost, const family_line_item& item)
{
  Report::Formats format = static_cast<Report::Formats>(ost.iword(Report::report_format));
  switch (format)
  {
  case Report::HTML:
    ost << item.node_id << "\t";
    ost << item.tree << "\t";
    ost << item.max_p_value << "\t(";
    for (size_t b = 0; b < item.pvalues.size(); b++)
    {
      std::pair<double, double> p = item.pvalues[b];
      if (p.first < 0)
      {
        ost << "(-,-)";
      }
      else
      {
        ost << p;
      }
      if (b < item.pvalues.size() - 1)
        ost << ",";
    }
    ost << ")\t";

    if (!item.cut_pvalues.empty())
    {
      write_doubles(ost, item.cut_pvalues);
      ost << "\t";
    }

    if (!item.likelihood_ratios.empty())
    {
      write_doubles(ost, item.likelihood_ratios);
    }
    break;
  case Report::JSON:
    ost << "{ " << quoted("Node") << ":" << quoted(item.node_id) << "," << endl;
    ost << quoted("Tree") << ":" << quoted(item.tree) << "," << endl;
    ost << quoted("MaxPValue") << ":" << item.max_p_value << "," << endl;
    ost << quoted("PValues") << ":[" << item.pvalues << "]}" << endl;
    break;
  case Report::Text:
  case Report::Unknown:
    ost << item.node_id << "\t";
    ost << item.tree << "\t";
    ost << item.max_p_value << "\t(";
    for (size_t b = 0; b < item.pvalues.size(); b++)
    {
      std::pair<double, double> p = item.pvalues[b];
      if (p.first < 0)
      {
        ost << "(-,-)";
      }
      else
      {
        ost << "(" << p.first << "," << p.second << ")";
      }
      if (b < item.pvalues.size() - 1)
        ost << ",";
    }
    ost << ")\t";

    if (!item.cut_pvalues.empty())
    {
      write_doubles(ost, item.cut_pvalues);
      ost << "\t";
    }

    if (!item.likelihood_ratios.empty())
    {
      write_doubles(ost, item.likelihood_ratios);
    }
    break;
  }

	return ost;
}

Report::Report(pCafeParam param, viterbi_parameters& viterbi)
{
	pString pstr = phylogeny_string((pTree)param->pcafe, NULL);
	tree = pstr->buf;
	string_free(pstr);

	for (int i = 0; i < param->num_lambdas; i++)
	{
		lambdas.push_back(param->lambda[i]);
	}

	if (param->lambda_tree)
	{
		pString pstr = phylogeny_string_newick(param->lambda_tree, lambda_tree_string, 0);
		lambda_tree = pstr->buf;
		string_free(pstr);
	}

	averageExpansion = viterbi.averageExpansion;
	changes = viterbi.expandRemainDecrease;

	pstr = phylogeny_string_newick((pTree)param->pcafe, cafe_tree_string_id, PS_SKIP_BL);
	id_tree = pstr->buf;
	string_free(pstr);

	pArrayList nlist = param->pcafe->super.nlist;

	for (int b = 1; b < nlist->size; b += 2)
	{
		pTreeNode pnode = (pTreeNode)nlist->array[b];
		pTreeNode child[2] = { (pTreeNode)tree_get_child(pnode,0), (pTreeNode)tree_get_child(pnode,1) };
		node_pairs.push_back(std::pair<int, int>(child[0]->id, child[1]->id));
	}

	pArrayList pflist = param->pfamily->flist;
	for (int i = 0; i < pflist->size; i++)
	{
		family_line_item item(param->pfamily, param->pcafe, param->likelihoodRatios, viterbi, i, ((pCafeFamilyItem)pflist->array[i])->id);
		family_line_items.push_back(item);
	}

	branch_cutting_output_format.resize(nlist->size);
	for (int i = 0; i < nlist->size; i++)
	{
		branch_cutting_output_format[i] = i;
	}
}

std::ostream& operator<<(ostream& ost, const Report& report)
{
  Report::Formats format = static_cast<Report::Formats>(ost.iword(Report::report_format));
  bool has_pvalues = !report.family_line_items.empty() && report.family_line_items[0].cut_pvalues.empty();
  bool has_likelihoods = !report.family_line_items.empty() && report.family_line_items[0].likelihood_ratios.empty();

  switch (format)
  {
  case Report::Text:
  case Report::Unknown:
    ost << "Tree:";
    ost << report.tree << "\n" << report.tree;

    ost << "Lambda:";
    for (size_t i = 0; i < report.lambdas.size(); i++)
    {
      ost << "\t" << report.lambdas[i];
    }
    ost << "\n";
    if (!report.lambda_tree.empty())
    {
      ost << "Lambda tree:\t" << report.lambda_tree << "\n";
    }

    ost << "# IDs of nodes:";

    ost << report.id_tree << "\n";

    ost << "# Output format for: ' Average Expansion', 'Expansions', 'No Change', 'Contractions', and 'Branch-specific P-values' = (node ID, node ID): ";
    for (size_t i = 0; i < report.node_pairs.size(); ++i)
    {
      ost << report.node_pairs[i] << " ";
    }

    ost << "\n";
    ost << "# Output format for 'Branch cutting P-values' and 'Likelihood Ratio Test': (" << report.branch_cutting_output_format[0];

    for (size_t i = 1; i < report.branch_cutting_output_format.size(); i++)
    {
      ost << ", " << report.branch_cutting_output_format[i];
    }
    ost << ")\n";

    write_viterbi(ost, report);

    write_families_header(ost, has_pvalues, has_likelihoods);

    copy(report.family_line_items.begin(), report.family_line_items.end(), ostream_iterator<family_line_item>(ost, "\n"));
    break;
  case Report::HTML:
    ost << "<!DOCTYPE html><html><head><meta charset = \"UTF-8\">";
    ost << "<title>Cafe Report</title><style>table, th, td{ border: 1px solid black;}</style></head><body>";
    ost << "<h1>Cafe Report</h1>";
    ost << "<h3>Input tree</h3><p>" << report.tree << "</p>";
    ost << "<h3>Lambda values</h3><ol>";
    for (size_t i = 0; i < report.lambdas.size(); i++)
    {
      ost << "<li>" << report.lambdas[i] << "</li>";
    }
    ost << "</ol>";
    if (!report.lambda_tree.empty())
    {
      ost << "<h3>Lambda tree</h3><p>" << report.lambda_tree << "</p>";
    }
    ost << "<h3>ID tree</h3><p>" << report.id_tree << "</p>";

    ost << "<h3>Family Change Summary</h3>";
    ost << "<table><tr><th>Node Pairs</th>";

    for (size_t i = 0; i < report.node_pairs.size(); ++i)
    {
      ost << "<th>" << report.node_pairs[i] << "</th>";
    }
    ost << "</tr><tr>";

    write_viterbi(ost, report);
    ost << "</tr></table>";

    ost << "<h3>Families</h3><table><tr><th>ID</th><th>Newick</th><th>Family-wide P-value</th><th>Viterbi P-values</th>";

    if (has_pvalues)
      ost << "<th>cut P-value</th>";
    if (has_likelihoods)
      ost << "<th>Likelihood Ratio</th>";
    ost << "</tr>";

    for (size_t i = 0; i < report.family_line_items.size(); ++i)
    {
      const family_line_item& item = report.family_line_items[i];
      ost << "<tr><td>" << item.node_id << "</td><td>" << item.tree << "</td><td>" << item.max_p_value << "</td><td>";
      for (size_t b = 0; b < item.pvalues.size(); b++)
      {
        std::pair<double, double> p = item.pvalues[b];
        if (p.first < 0)
        {
          ost << "(-,-)";
        }
        else
        {
          ost << "(" << p.first << "," << p.second << ")";
        }
        if (b < item.pvalues.size() - 1)
          ost << ",";
      }
      ost << "</td>";

      if (!item.cut_pvalues.empty())
      {
        ost << "<td>";
        write_doubles(ost, item.cut_pvalues);
        ost << "</td>";
      }

      if (!item.likelihood_ratios.empty())
      {
        ost << "<td>";
        write_doubles(ost, item.likelihood_ratios);
        ost << "</td>";
      }
      ost << "</tr>";
    }
    ost << "</table>";

    ost << "</body></html>";
    break;
  case Report::JSON:
    ost << json;

    ost << "{";
    ost << quoted("InputTree") << ":" << quoted(report.tree) << "," << endl;
    ost << quoted("Lambdas") << ":[" << report.lambdas << "]," << endl;
    if (!report.lambda_tree.empty())
    {
      ost << quoted("LambdaTree") << ":" << quoted(report.lambda_tree) << "," << endl;
    }
    ost << quoted("IDTree") << ":" << quoted(report.id_tree) << "," << endl;
    ost << quoted("NodePairs") << ":[" << report.node_pairs << "]," << endl;
    write_viterbi(ost, report);
    ost << quoted("Families") << ":[" << report.family_line_items << "]" << endl;
    ost << "}";
    break;
  }

	return ost;
}

report_parameters get_report_parameters(std::vector<std::string> tokens)
{
	report_parameters params;
	params.name = tokens[1];

	params.branchcutting = false;
	params.likelihood = false;
	params.lh2 = false;
	params.just_save = false;
  params.format = Report::Text;

	for (size_t i = 2; i < tokens.size(); i++)
	{
		if (strcasecmp(tokens[i].c_str(), "html") == 0) params.format = Report::HTML;
    if (strcasecmp(tokens[i].c_str(), "json") == 0) params.format = Report::JSON;
    if (strcasecmp(tokens[i].c_str(), "branchcutting") == 0) params.branchcutting = true;
		if (strcasecmp(tokens[i].c_str(), "likelihood") == 0) params.likelihood = true;
		if (strcasecmp(tokens[i].c_str(), "lh2") == 0) params.lh2 = true;
		if (strcasecmp(tokens[i].c_str(), "save") == 0)
		{
			params.branchcutting = false;
			params.likelihood = false;
			params.lh2 = false;
			params.just_save = true;
			break;
		}
	}
	return params;

}

string extension(Report::Formats f)
{
  if (f == Report::HTML)
    return ".html";
  if (f == Report::JSON)
    return ".json";
  if (f == Report::Text)
    return ".cafe";

  throw runtime_error("Unknown Report Type");
}

void set_format(ostream& ost, Report::Formats f)
{
  if (f == Report::HTML)
    ost << html;
  else if (f == Report::JSON)
    ost << json;
}

void cafe_do_report(Globals& globals, viterbi_parameters& viterbi, report_parameters* params)
{
	pCafeParam param = &globals.param;
	if (!params->just_save)
	{
		cafe_tree_set_parameters(param->pcafe, &param->family_size, 0);
		int nnodes = ((pTree)param->pcafe)->nlist->size;
    viterbi.clear(nnodes);
	}

  string filename = params->name + extension(params->format);
	ofstream report(filename.c_str());
  set_format(report, params->format);

	if (!report)
	{
		throw std::runtime_error(string("ERROR(report) : Cannot open ") + params->name + " in write mode.\n");
	}

	if (ConditionalDistribution::matrix.empty())
	{
		param->param_set_func(param, param->input.parameters);
		reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, &param->family_size);
		ConditionalDistribution::reset(param->pcafe, &param->family_size, param->num_threads, globals.num_random_samples);
	}


	if (params->branchcutting || params->likelihood)
	{
		pArrayList cd = ConditionalDistribution::to_arraylist();
		cafe_viterbi(globals, viterbi, cd);
		arraylist_free(cd, NULL);
		
		if (params->branchcutting) 
			cafe_branch_cutting(globals, viterbi);
		
		if (params->likelihood) 
			cafe_likelihood_ratio_test(param, viterbi.maximumPvalues);
		
		cafe_log(param, "Building Text report: %s\n", params->name.c_str());
		Report r(param, viterbi);
		report << r;
	}
	else if (params->lh2)
	{
		cafe_lhr_for_diff_lambdas(param, tmp_lambda_tree, 2, cafe_shell_set_lambda);
	}
	else
	{
		if (!params->just_save)
		{
			pArrayList cd = ConditionalDistribution::to_arraylist();
			cafe_viterbi(globals, viterbi, cd);
			arraylist_free(cd, NULL);
		}
		Report r(param, viterbi);
		report << r;
	}


	cafe_log(param, "Report Done\n");

}

void cafe_report_load_viterbi_pvalue(char* data, viterbi_parameters& v, pCafeFamilyItem item, pCafeTree ptree)
{
  int nnodes = ptree->super.nlist->size;

  char* next = &data[1];
  while ((next = strchr(next, ')')) != NULL)
  {
    next++;
    if (next[0] == ',')
    {
      next[0] = '\t';
    }
    else if (next[0] == ')')
    {
      next[0] = '\0';
      break;
    }
  }
  double* array = cafe_report_load_data_double_pairs(&data[1], '\t');
  int j;
  for (j = 0; j < nnodes; j++)
  {
    pCafeNode pnode = (pCafeNode)ptree->super.nlist->array[j];
    v.viterbiPvalues[viterbi_parameters::NodeFamilyKey(pnode, item)] = array[j];
  }
  memory_free(array);
  array = NULL;
}


int cafe_report_retrieve_data(const char* file, pCafeParam param, viterbi_parameters& viterbi)
{
	int i, j;
	family_size_range range;
	range.min = 0;
	range.max = 1;
	range.root_min = 1;
	range.root_max = 2;

	FILE* fp = fopen(file, "r");
	if (fp == NULL)
	{
		fprintf(stderr, "%s: Cannot open file: %s\n", __FUNCTION__, file);
		return -1;
	}
	pString pstr = string_new();
	int nnodes = 0;
	viterbi.expandRemainDecrease.clear();

	int bexist[2] = { 0 , 0 };

	while (file_read_line(pstr, fp))
	{
		if (pstr->buf[0] == '\'')
		{
			pArrayList list = string_pchar_split(pstr->buf, '\t');
			if (list->size == 4) break;
			else if (list->size == 6)
			{
				bexist[0] = 1;
				bexist[1] = 1;
			}
			else
			{
				if ((pstr->buf[1] & 0x4F) == 'C')
				{
					bexist[0] = 1;
				}
				else
				{
					bexist[1] = 1;
				}
			}
			break;
		}
		char* data = index(pstr->buf, ':');
		if (data == NULL)
		{
			fprintf(stderr, "Format error during loading cafe file: %s\n", pstr->buf);
			fclose(fp);
			return -2;
		}
		*data++ = '\0';
		string_pchar_chomp(data);
		if (strcasecmp(pstr->buf, "tree") == 0)
		{
			param->pcafe = cafe_tree_new(data, &range, 0, 0);
			nnodes = param->pcafe->super.nlist->size;
		}
		else if (strncasecmp(pstr->buf, "lambda tree", 10) == 0)
		{
			param->lambda_tree = phylogeny_new(data, phylogeny_lambda_parse_func);
			tree_build_node_list(param->lambda_tree);
		}
		else if (strncasecmp(pstr->buf, "lambda", 5) == 0)
		{
			param->lambda = cafe_report_load_data_double_list(data, '\t', &param->num_lambdas);
		}
		else if (strncasecmp(pstr->buf, "average", 7) == 0)
		{
			pArrayList list = string_pchar_split(data, '\t');
			viterbi.averageExpansion.resize(list->size * 2);
			int i;
			for (i = 0; i < list->size; i++)
			{
				viterbi.averageExpansion[2 * i] = -1;
				viterbi.averageExpansion[2 * i + 1] = -1;
				sscanf((char*)list->array[i], "(%lf,%lf)", &viterbi.averageExpansion[i * 2], &viterbi.averageExpansion[i * 2 + 1]);
			}
			arraylist_free(list, free);
		}
		else if (strncasecmp(pstr->buf, "expansion", 9) == 0)
		{
			pArrayList list = string_pchar_split(data, '\t');
			int sz = max(viterbi.expandRemainDecrease.size(), (size_t)(list->size + 1) * 2);
			viterbi.expandRemainDecrease.resize(sz);
			for (int i = 0; i < list->size; i++)
			{
				viterbi.expandRemainDecrease[i * 2].expand = -1;
				viterbi.expandRemainDecrease[i * 2 + 1].expand = -1;
				sscanf((char*)list->array[i], "(%d,%d)", &viterbi.expandRemainDecrease[i * 2].expand, &viterbi.expandRemainDecrease[i * 2 + 1].expand);
			}
			arraylist_free(list, free);
		}
		else if (strncasecmp(pstr->buf, "remain", 6) == 0)
		{
			pArrayList list = string_pchar_split(data, '\t');
			int sz = max(viterbi.expandRemainDecrease.size(), (size_t)(list->size + 1) * 2);
			viterbi.expandRemainDecrease.resize(sz);
			int i;
			for (i = 0; i < list->size; i++)
			{
				viterbi.expandRemainDecrease[i * 2].remain = -1;
				viterbi.expandRemainDecrease[i * 2 + 1].remain = -1;
				sscanf((char*)list->array[i], "(%d,%d)", &viterbi.expandRemainDecrease[i * 2].remain, &viterbi.expandRemainDecrease[i * 2 + 1].remain);
			}
			arraylist_free(list, free);
		}
		else if (strncasecmp(pstr->buf, "decrease", 8) == 0)
		{
			pArrayList list = string_pchar_split(data, '\t');
			int sz = max(viterbi.expandRemainDecrease.size(), (size_t)(list->size + 1) * 2);
			viterbi.expandRemainDecrease.resize(sz);
			int i;
			for (i = 0; i < list->size; i++)
			{
				viterbi.expandRemainDecrease[i * 2].decrease = -1;
				viterbi.expandRemainDecrease[i * 2 + 1].decrease = -1;
				sscanf((char*)list->array[i], "(%d,%d)", &viterbi.expandRemainDecrease[i * 2].decrease, &viterbi.expandRemainDecrease[i * 2 + 1].decrease);
			}
			arraylist_free(list, free);
		}
	}
	param->param_set_func(param, param->lambda);
	pArrayList plines = arraylist_new(11000);
	int num_families;
	for (num_families = 0; file_read_line(pstr, fp); num_families++)
	{
		char* line = (char*)memory_new(pstr->length + 1, sizeof(char));
		strcpy(line, pstr->buf);
		arraylist_add(plines, line);
	}

	viterbi_parameters_init(&viterbi, nnodes, num_families);

	param->pfamily = (pCafeFamily)memory_new(1, sizeof(CafeFamily));
	pCafeFamily pcf = param->pfamily;
	pcf->flist = arraylist_new(11000);
	pcf->num_species = (nnodes + 1) / 2;
	pcf->species = (char**)memory_new(pcf->num_species, sizeof(char*));
	pcf->index = (int*)memory_new(pcf->num_species, sizeof(int));
	pArrayList nlist = param->pcafe->super.nlist;


	for (i = j = 0; i < nnodes; i += 2, j++)
	{
		char* name = ((pPhylogenyNode)nlist->array[i])->name;
		pcf->index[j] = i;
		pcf->species[j] = (char*)memory_new(strlen(name) + 1, sizeof(char));
		strcpy(pcf->species[j], name);
	}

	if (bexist[0])
	{
		viterbi.cutPvalues = (double**)memory_new_2dim(nnodes, num_families, sizeof(double));
	}
	if (bexist[1])
	{
		param->likelihoodRatios = (double**)memory_new_2dim(nnodes, num_families, sizeof(double));
	}
	int max_size = 0;
	for (i = 0; i < plines->size; i++)
	{
		pArrayList data = string_pchar_split((char*)plines->array[i], '\t');
		pCafeFamilyItem pitem = (pCafeFamilyItem)memory_new(1, sizeof(CafeFamilyItem));
		pitem->desc = NULL;
		pitem->id = (char*)data->array[0];
		pitem->count = (int*)memory_new(pcf->num_species, sizeof(int));
		pitem->maxlh = -1;
		pitem->ref = -1;
		pitem->lambda = NULL;
		arraylist_add(pcf->flist, pitem);
		sscanf((char*)data->array[2], "%lf", &viterbi.maximumPvalues[i]);

		pCafeTree ptree = cafe_tree_new((char*)data->array[1], &range, 0, 0);
		pArrayList nlist = ptree->super.nlist;

		for (j = 0; j < nlist->size; j += 2)
		{
			pCafeNode pcnode = (pCafeNode)nlist->array[j];
			pitem->count[j / 2] = pcnode->familysize;
			if (max_size < pcnode->familysize)
			{
				max_size = pcnode->familysize;
			}
		}
		for (j = 1; j < nlist->size; j += 2)
		{
			pCafeNode pcnode = (pCafeNode)nlist->array[j];
			viterbi.viterbiNodeFamilysizes[viterbi_parameters::NodeFamilyKey(pcnode, pitem)] = pcnode->familysize;
		}

		cafe_report_load_viterbi_pvalue((char*)data->array[3], viterbi, pitem, ptree);
		if (bexist[0])
		{
			cafe_report_load_bc_or_lhr_list((char*)data->array[4],
				viterbi.cutPvalues, i, nnodes);
			viterbi.cutPvalues[param->pcafe->super.root->id][i] = -1.0;
		}
		if (bexist[1])
		{
			int idx = bexist[0] ? 5 : 4;
			cafe_report_load_bc_or_lhr_list((char*)data->array[idx],
				param->likelihoodRatios, i, nnodes);
		}
		cafe_tree_free(ptree);
	}

	init_family_size(&range, max_size);

	cafe_tree_set_parameters(param->pcafe, &range, param->lambda[0]);
	arraylist_free(plines, free);
	return 0;
}



