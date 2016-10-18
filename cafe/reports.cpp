#include <stdexcept>
#include "reports.h"

using namespace std;

extern "C" {
#include "viterbi.h"
#include "cafe.h"

	extern pTree tmp_lambda_tree;
	void cafe_shell_set_lambda(pCafeParam param, double* lambda);
}

pArrayList cafe_pCD;

void cafe_do_report(pCafeParam param, report_parameters* params)
{
	if (!params->just_save)
	{
		cafe_shell_set_sizes();
		int nnodes = ((pTree)param->pcafe)->nlist->size;
		viterbi_parameters_clear(&param->viterbi, nnodes);
	}

	if (params->bc || params->lh)
	{
		cafe_pCD = cafe_viterbi(param, cafe_pCD);
		if (params->bc) cafe_branch_cutting(param);
		if (params->lh) cafe_likelihood_ratio_test(param);
		cafe_log(param, "Building Text report: %s\n", params->name.c_str());
		cafe_report(param, CAFE_REPORT_TEXT);
		fclose(param->fout);
		string name = params->name + ".mp";
		param->fout = fopen(name.c_str(), "w");
		cafe_log(param, "Building Metapost report: %s\n", name.c_str());
		cafe_report(param, CAFE_REPORT_PDF);
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
		cafe_report(param, CAFE_REPORT_TEXT);
		fclose(param->fout);
		string name = params->name + ".mp";
		param->fout = fopen(name.c_str(), "w");
		cafe_log(param, "Building Metapost report: %s\n", name.c_str());
		cafe_report(param, CAFE_REPORT_PDF);
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

int cafe_cmd_report(pCafeParam param, std::vector<std::string> tokens)
{
	if (param->pfamily == NULL)
		throw std::runtime_error("ERROR(report): You did not load family: command 'load'\n");
	if (param->pcafe == NULL)
		throw std::runtime_error("ERROR(report): You did not specify tree: command 'tree'\n");
	if (param->lambda == NULL)
		throw std::runtime_error("ERROR(report): You did not set the parameters: command 'lambda' or 'lambdamu'\n");

	report_parameters params;
	get_report_parameters(params, tokens);

	string name = params.name + ".cafe";
	param->fout = fopen(name.c_str(), "w");

	if (param->fout == NULL)
	{
		throw std::runtime_error(string("ERROR(report) : Cannot open ") + name + " in write mode.\n");
	}

	cafe_do_report(param, &params);
	return 0;
}



