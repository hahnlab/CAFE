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
 
void cafe_shell_prompt(char* prompt, char* format, ... );
double cafe_shell_score();
int cafe_cmd_lambda_mu(int argc, char* argv[]);
int cafe_cmd_esterror(int argc, char* argv[]);
int cafe_cmd_simerror(int argc, char* argv[]);

int set_log_file(pCafeParam param, const char *file_name);

void initialize_k_bd(pCafeParam param, double *parameters);
void set_birth_death_probabilities4(struct probabilities *probs, int num_lambdas, int fix_cluster, int taxa_id, double* parameters);
void initialize_k_weights(pCafeParam param);
void initialize_k_weights2(pCafeParam param);

#endif 
