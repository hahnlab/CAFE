#ifndef VITERBI_H_A08989A1_B4B4_461C_B863_A1AE2FE9BD98
#define VITERBI_H_A08989A1_B4B4_461C_B863_A1AE2FE9BD98

struct viterbi_parameters;

void viterbi_set_max_pvalue(viterbi_parameters* viterbi, int index, double val);
void viterbi_parameters_init(viterbi_parameters *viterbi, int nnodes, int nrows);
void viterbi_parameters_clear(viterbi_parameters* viterbi, int nnodes);

#endif
