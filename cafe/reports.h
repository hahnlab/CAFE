#ifndef REPORTS_C6B4B223_5C6D_459E_B2D8_CAC46A40A5B4
#define REPORTS_C6B4B223_5C6D_459E_B2D8_CAC46A40A5B4

#include <vector>
#include <string>

extern "C" {
#include "cafe_shell.h"
#include "family.h"
}

typedef struct
{
	int bc, lh, lh2, just_save;
	std::string name;
} report_parameters;


void get_report_parameters(report_parameters &params, std::vector<std::string> tokens);
int cafe_cmd_report(pCafeParam param, std::vector<std::string> tokens);

#endif
