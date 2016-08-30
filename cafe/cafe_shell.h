#ifndef __CAFE_SHELL_H__
#define __CAFE_SHELL_H__

#include<stdarg.h>
#include<stdio.h>
#include<cafe.h>
#include<utils.h>
#include<io.h>

typedef int (*cafe_command)(int argc, char* argv[]);

typedef struct 
{
	char* command; 
	cafe_command func;
}CafeShellCommand;

typedef CafeShellCommand* pCafeShellCommand;
extern CafeShellCommand cafe_cmd[];
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
int cafe_shell_dispatch_command(char* cmd);
void cafe_shell_prompt(char* prompt, char* format, ... );
void cafe_shell_init();
double cafe_shell_score();
int cafe_cmd_tree(int argc, char* argv[]);
int cafe_cmd_branchlength(int argc, char* argv[]);
int cafe_cmd_echo(int argc, char* argv[]);
int cafe_cmd_family(int argc, char* argv[]);
int cafe_cmd_list(int argc, char* argv[]);
int cafe_cmd_date(int argc, char* argv[]);
int cafe_cmd_viterbi(int argc, char* argv[]);
int cafe_cmd_print_param(int argc, char* argv[]);
int cafe_cmd_load(int argc, char* argv[]);
int cafe_cmd_save(int argc, char* argv[]);
int cafe_cmd_lambda(int argc, char* argv[]);
int cafe_cmd_lambda_mu(int argc, char* argv[]);
int cafe_cmd_version(int argc, char* argv[]);
int cafe_cmd_pvalue(int argc, char* argv[]);
int cafe_cmd_score(int argc, char* argv[]);
int cafe_cmd_report(int argc, char* argv[] );
int cafe_cmd_log(int argc, char* argv[] );
int cafe_cmd_cgi(int argc, char* argv[] );
int cafe_cmd_retrieve(int argc, char* argv[] );
int cafe_cmd_gainloss(int argc, char* argv[] );
int cafe_cmd_sim_extinct(int argc, char* argv[]);
int cafe_cmd_extinct(int argc, char* argv[]);
int cafe_cmd_generate_random_family(int argc, char* argv[]);
int cafe_cmd_lh_test(int argc, char* argv[]);
int cafe_cmd_root_dist(int argc, char* argv[]);
int cafe_cmd_esterror(int argc, char* argv[]);
int cafe_cmd_error_model(int argc, char* argv[]);
int cafe_cmd_no_error_model(int argc, char* argv[]);
int cafe_cmd_simerror(int argc, char* argv[]);
int cafe_cmd_reconstruction_accuracy(int argc, char* argv[]);
int cafe_cmd_crossvalidation_by_species(int argc, char* argv[]);
int cafe_cmd_crossvalidation_by_family(int argc, char* argv[]);
int cafe_cmd_source(int argc, char* argv[]);
int cafe_cmd_exit(int argc, char* argv[]);



#endif 
