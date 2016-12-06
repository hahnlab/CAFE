#include <stdexcept>
#include <ostream>
#include <iomanip>
#include <fstream>

#include "reports.h"
#include "likelihood_ratio.h"

using namespace std;

extern "C" {
#include "viterbi.h"
#include "cafe.h"

	extern pTree tmp_lambda_tree;
	void cafe_shell_set_lambda(pCafeParam param, double* lambda);
}

pArrayList cafe_pCD;

void lambda_tree_string(pString pstr, pPhylogenyNode pnode)
{
	if (pnode->taxaid != -1)
	{
		string_fadd(pstr, "%d", pnode->taxaid + 1);
	}
}

void cafe_report_set_viterbi(pCafeParam param, int i)
{
	int j;
	cafe_family_set_size(param->pfamily, i, param->pcafe);
	pArrayList nlist = param->pcafe->super.nlist;
	for (j = 1; j < nlist->size; j += 2)
	{
		pCafeNode pcnode = (pCafeNode)nlist->array[j];
		pcnode->familysize = param->viterbi.viterbiNodeFamilysizes[j / 2][i];
	}
}

void cafe_tree_string_id(pString pstr, pPhylogenyNode pnode)
{
	if (pnode->name) string_fadd(pstr, "%s", pnode->name);
	string_fadd(pstr, "<%d>", pnode->super.id);
}


void write_viterbi(ostream& ost, const viterbi_parameters& viterbi)
{
	ost << "Average Expansion:";
	int nnodes = viterbi.num_nodes / 2;
	for (int b = 0; b < nnodes; b++)
	{
		ost << "\t(" << viterbi.averageExpansion[2 * b] << "," << viterbi.averageExpansion[2 * b + 1] << ")";
	}
	ost << "\n";
	ost << "Expansion :";
	for (int b = 0; b < nnodes; b++)
	{
		ost << "\t(" << viterbi.expandRemainDecrease[0][2 * b] << "," << viterbi.expandRemainDecrease[0][2 * b + 1] << ")";
	}
	ost << "\n";

	ost << "nRemain :";
	for (int b = 0; b < nnodes; b++)
	{
		ost << "\t(" << viterbi.expandRemainDecrease[1][2 * b] << "," << viterbi.expandRemainDecrease[1][2 * b + 1] << ")";
	}
	ost << "\n";

	ost << "nDecrease :";
	for (int b = 0; b < nnodes; b++)
	{
		ost << "\t(" << viterbi.expandRemainDecrease[2][2 * b] << "," << viterbi.expandRemainDecrease[2][2 * b + 1] << ")";
	}
	ost << "\n";
}

void write_families_header(ostream& ost, double **cutPvalues, double**likelihoodRatios)
{
	ost << "'ID'\t'Newick'\t'Family-wide P-value'\t'Viterbi P-values'";
	if (cutPvalues)
		ost << "\t'cut P-value'";
	if (likelihoodRatios)
		ost << "\t'Likelihood Ratio'";
	ost << "\n";
}

void write_families_line(ostream& ost, pCafeParam param, int i, string node_id)
{
	ost << node_id << "\t";
	cafe_report_set_viterbi(param, i);
	pString pstr = cafe_tree_string(param->pcafe);
	ost << pstr->buf << "\t";
	string_free(pstr);
	ost << param->viterbi.maximumPvalues[i] << "\t(";
	for (int b = 0; b < param->viterbi.num_nodes / 2; b++)
	{
		if (param->viterbi.viterbiPvalues[2 * b][i] == -1)
		{
			ost << "(-,-)";
		}
		else
		{
			ost << "(" << param->viterbi.viterbiPvalues[2 * b][i] << "," << param->viterbi.viterbiPvalues[2 * b + 1][i] << ")";
		}
		if (b < param->viterbi.num_nodes / 2 - 1)
			ost << ",";
	}
	ost << ")\t";

	if (param->viterbi.cutPvalues)
	{
		ost << "(";
		for (int b = 0; b < param->viterbi.num_nodes; b++)
		{
			if (param->viterbi.cutPvalues[b][i] == -1)
			{
				ost << "-";
			}
			else
			{
				ost << param->viterbi.cutPvalues[b][i];
			}
			if (b < param->viterbi.num_nodes - 1) 
				ost << ",";
		}
		ost << ")\t";
	}

	if (param->likelihoodRatios)
	{
		ost << "(";
		for (int b = 0; b < param->pcafe->super.nlist->size; b++)
		{
			if (param->likelihoodRatios[b][i] == -1)
			{
				ost << "-";
			}
			else
			{
				ost << param->likelihoodRatios[b][i];
			}
			if (b <param->pcafe->super.nlist->size - 1) 
				ost << ",";
		}
		ost << ")";
	}

	ost << "\n";

}

void cafe_report(pCafeParam param, ostream& report_file)
{
	report_file << "Tree:";
	pString pstr = phylogeny_string((pTree)param->pcafe, NULL);
	report_file << pstr->buf << "\n" << pstr->buf;
	string_free(pstr);

	report_file << "Lambda:";
	for (int i = 0; i < param->num_lambdas; i++)
	{
		report_file << "\t" << param->lambda[i];
	}
	report_file << "\n";
	if (param->lambda_tree)
	{
		pString pstr = phylogeny_string_newick(param->lambda_tree, lambda_tree_string, 0);
		report_file << "Lambda tree:\t" << pstr->buf << "\n";
		string_free(pstr);
	}

	report_file << "# IDs of nodes:";
	//	pstr = cafe_tree_string_with_id(param->pcafe);

	pstr = phylogeny_string_newick((pTree)param->pcafe, cafe_tree_string_id, PS_SKIP_BL);
	report_file << pstr->buf << "\n";
	string_free(pstr);

	pArrayList nlist = param->pcafe->super.nlist;

	report_file << "# Output format for: ' Average Expansion', 'Expansions', 'No Change', 'Contractions', and 'Branch-specific P-values' = (node ID, node ID): ";
	for (int b = 1; b < nlist->size; b += 2)
	{
		pTreeNode pnode = (pTreeNode)nlist->array[b];
		pTreeNode child[2] = { (pTreeNode)tree_get_child(pnode,0), (pTreeNode)tree_get_child(pnode,1) };
		report_file << "(" << child[0]->id << "," << child[1]->id << ") ";
	}

	report_file << "\n";
	report_file << "# Output format for 'Branch cutting P-values' and 'Likelihood Ratio Test': (" << 0;
	for (int i = 1; i < nlist->size; i++)
	{
		report_file << ", " << i;
	}
	report_file << ")\n";

	write_viterbi(report_file, param->viterbi);

	write_families_header(report_file, param->viterbi.cutPvalues, param->likelihoodRatios);

	pArrayList pflist = param->pfamily->flist;
	for (int i = 0; i < pflist->size; i++)
	{
		write_families_line(report_file, param, i, ((pCafeFamilyItem)pflist->array[i])->id);
	}
}

void cafe_do_report(pCafeParam param, report_parameters* params)
{
	if (!params->just_save)
	{
		cafe_shell_set_sizes();
		int nnodes = ((pTree)param->pcafe)->nlist->size;
		viterbi_parameters_clear(&param->viterbi, nnodes);
	}

	string filename = params->name + ".cafe";
	ofstream report(filename.c_str());

	if (!report)
	{
		throw std::runtime_error(string("ERROR(report) : Cannot open ") + params->name + " in write mode.\n");
	}

	if (params->bc || params->lh)
	{
		cafe_pCD = cafe_viterbi(param, cafe_pCD);
		if (params->bc) cafe_branch_cutting(param);
		if (params->lh) cafe_likelihood_ratio_test(param);
		cafe_log(param, "Building Text report: %s\n", params->name.c_str());
		cafe_report(param, report);
	}
	else if (params->lh2)
	{
		cafe_lhr_for_diff_lambdas(param, tmp_lambda_tree, 2, cafe_shell_set_lambda);
	}
	else
	{
		if (!params->just_save)
		{
			cafe_pCD = cafe_viterbi(param, cafe_pCD);
		}
		cafe_report(param, report);
	}
	fclose(param->fout);

	// HTML
	string name = params->name + ".html";
	FILE* fhttp = fopen(name.c_str(), "w");

	cafe_log(param, "Building HTML report: %s\n", name.c_str());
	fprintf(fhttp, "<html>\n<body>\n<table border=1>\n");
	for (int i = 0; i < param->pfamily->flist->size; i++)
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
		fprintf(fhttp, "<tr><td><a href=pdf/%s-%d.pdf>%s</a></td><td>%s</td></tr>\n",
			params->name.c_str(), i + 1, pitem->id, pitem->desc ? pitem->desc : "NONE");
	}
	fprintf(fhttp, "</table>\n</body>\n</html>\n");
	fclose(fhttp);

	cafe_log(param, "Report Done\n");

}

void get_report_parameters(report_parameters &params, std::vector<std::string> tokens)
{
	params.name = tokens[1];

	params.bc = 0;
	params.lh = 0;
	params.lh2 = 0;
	params.just_save = 0;
	for (size_t i = 2; i < tokens.size(); i++)
	{
		if (strcasecmp(tokens[i].c_str(), "branchcutting") == 0) params.bc = 1;
		if (strcasecmp(tokens[i].c_str(), "likelihood") == 0) params.lh = 1;
		if (strcasecmp(tokens[i].c_str(), "lh2") == 0) params.lh2 = 1;
		if (strcasecmp(tokens[i].c_str(), "save") == 0)
		{
			params.bc = 0;
			params.lh = 0;
			params.lh2 = 0;
			params.just_save = 1;
			break;
		}
	}
}




