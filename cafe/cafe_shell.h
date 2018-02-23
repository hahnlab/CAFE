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
 
void initialize_k_bd(pCafeTree pcafe, pTree lambda_tree, int num_values, int fixcluster, double *parameters);
void set_birth_death_probabilities4(struct probabilities *probs, int num_lambdas, int fix_cluster, int taxa_id, double* parameters);

#endif 
