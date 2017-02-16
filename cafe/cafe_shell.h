#ifndef __CAFE_SHELL_H__
#define __CAFE_SHELL_H__

#include <utils.h>
#include <utils_string.h>
#include <family.h>

typedef enum
{
	CAFE_SHELL_EXIT = 10000,	
	CAFE_SHELL_NO_COMMAND,
}enumCafeSehll;

typedef struct
{
	char* opt;
	int   argc; 
	char** argv;
}Argument;
typedef Argument* pArgument;

typedef struct
{
    double* sizeDist;
    int maxFamilySize;
    int** pairs;
    int model_parameter_number;
    int model_parameter_diff;
    int b_symmetric;
    int b_peakzero;
    double* estimates; // estimates index goes in order row (-diff) ... -1 0 1 ... (+diff) for asymmetric models and 0 1 ... (diff) for symmetric models. 

} ErrorMeasure;
typedef ErrorMeasure* pErrorMeasure;
 
pArrayList cafe_shell_build_argument(int argc, char* argv[]);
void cafe_shell_prompt(char* prompt, char* format, ... );
void cafe_shell_init(int quiet);
double cafe_shell_score();
int cafe_cmd_lambda_mu(int argc, char* argv[]);
int cafe_cmd_esterror(int argc, char* argv[]);
int cafe_cmd_simerror(int argc, char* argv[]);

int set_log_file(pCafeParam param, const char *file_name);
void log_param_values(pCafeParam param);

#endif 
