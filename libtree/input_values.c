#include <stddef.h>
#include <stdlib.h>
#include "input_values.h"

void input_values_init(input_values* vals)
{
#if 0
	vals->num_lambdas = 0;
	vals->lambdas = NULL;

	vals->num_k_weights = 0;
	vals->k_weights = NULL;
#else
	vals->parameters = NULL;
	vals->num_lambdas = 0;
	vals->num_k_weights = 0;
	vals->first_k_weight = 0;
	vals->num_mus = 0;
	vals->num_allocated = 0;
#endif
}

void input_values_construct(input_values* vals, int count)
{
	free(vals->parameters);
	vals->parameters = (double*)calloc(count, sizeof(double));
	vals->num_allocated = count;
}

void input_values_destruct(input_values* vals)
{
	free(vals->parameters);
	vals->parameters = NULL;
	vals->num_lambdas = 0;
	vals->num_k_weights = 0;
	vals->first_k_weight = 0;
	vals->num_mus = 0;
	vals->first_mu = 0;
	vals->num_allocated = 0;
}

void input_values_copy(input_values *dest, input_values* src)
{
#if 0
	if (dest->lambdas)
		free(dest->lambdas);

	if (dest->k_weights)
		free(dest->k_weights);

	dest->lambdas = calloc(src->num_lambdas, sizeof(double));
	dest->k_weights = calloc(src->num_k_weights, sizeof(double));

	dest->num_lambdas = src->num_lambdas;
	dest->num_k_weights = src->num_k_weights;
#else
	input_values_destruct(dest);
#endif
}

void input_values_set_lambdas(input_values *dest, double *src, int count)
{
	dest->num_lambdas = count;
	for (int i = 0; i < count; ++i)
		dest->parameters[i] = src[i];
}

void input_values_set_k_weights(input_values *dest, double *src, int start, int count)
{
	dest->num_k_weights = count;
	dest->first_k_weight = start;
	for (int i = 0; i < count; ++i)
		dest->parameters[start+i] = src[i];
}

void input_values_set_mus(input_values *dest, double *src, int start, int count)
{
	dest->num_mus = count;
	dest->first_mu = start;
	for (int i = 0; i < count; ++i)
		dest->parameters[start + i] = src[i];
}

void input_values_copy_weights(double *out, input_values *src, int start, int count)
{
	int i;
	double sumofweights = 0;
	for (i = 0; i < count - 1; i++) {
		out[i] = src->parameters[start + i];
		sumofweights += src->parameters[start + i];
	}
	out[i] = 1 - sumofweights;

}
