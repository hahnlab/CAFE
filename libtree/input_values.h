
typedef struct {
	int num_lambdas;

	int num_k_weights;
	int first_k_weight;

	int num_mus;
	int first_mu;
#if 0

	double *lambdas;

	double *k_weights;
#else
	double *parameters;
	int num_allocated;
#endif
} input_values;

void input_values_init(input_values* vals);
void input_values_copy(input_values *dest, input_values* src);
void input_values_construct(input_values* vals, int count);
void input_values_destruct(input_values* vals);
void input_values_set_lambdas(input_values *dest, double *src, int count);
void input_values_set_mus(input_values *dest, double *src, int start, int count);
void input_values_set_k_weights(input_values *dest, double *src, int start, int count);
void input_values_copy_weights(double *out, input_values *src, int start, int count);
