#ifndef __MATHFUNC_H__
#define __MATHFUNC_H__

#include<math.h>
#include<stdarg.h>
#include<stdio.h>

#ifdef	__cplusplus
extern "C" {
#endif

#define __EPS__		1e-8

#ifndef MAX
#define MAX(a,b)  ( (a) > (b) ? (a) : (b) )
#endif

#ifndef MIN
#define MIN(a,b)  ( (a) < (b) ? (a) : (b) )
#endif

#define MAX_DOUBLE 	1.7976931348623157e+308
#define MIN_DOUBLE  4.9e-300


typedef double (*math_func)(double* x, void* args);

typedef struct
{
	int maxiters;		
	int bymax;
	double rho, chi, psi, sigma;
	double tolx, tolf;
	double delta, zero_delta;

	int 	N, N1;		
	int 	iters;
	double** v;
	double* fv;
	double** vsort;
	double* x_mean;
	double* x_r;
	double* x_tmp;
	int*    idx;

	void* args;
	math_func eq;	
}FMinSearch;

typedef FMinSearch* pFMinSearch;

extern int __maxidx(double* data, int size );
extern double __max(double* data, int size );
extern double __min(double* data, int size );
extern double ipow(double val, int expo);

extern double unifrnd();

extern double gamma(double c);
extern double gammaln(double c);
extern double gammaincln(double x, double a);
extern double gammainc(double x, double a);

extern double beta(double alpha, double beta);
extern double betaln(double alpha, double beta);
extern double betainc(double x, double alpha, double beta);

extern double choose(double n, double r);
extern double chooseln(double n, double r);
extern double gampdf(double x, double alpha, double beta);
extern double gamcdf(double x, double alpha, double beta);
extern double betapdf(double x, double alpha, double beta);
extern double chi2pdf(double x, int df);
extern double chi2cdf(double x, int df);
extern double fpdf(double x, int df1, int df2 );
extern double fcdf(double x, int df1, int df2 );
extern double tpdf(double x, int df);
extern double tcdf(double x, int df);
extern double binopdf(int x, int n, double p);
extern double binocdf(int x, int n, double p);
extern double poisspdf(int x, double lambda);
extern double poisscdf(int x, double lambda );
extern double exppdf(double x, double mu );
extern double exprnd(double mu);
extern double reylpdf(double x, double b);
extern double normpdf(double x, double mu, double sigma);
extern double normrnd(double mu, double sigma);
extern double normcdf(double x, double mu, double sigma);
extern double norminv(double x, double mu, double sigma);
extern double mean(double* data, int size);
extern double variance(double* data, int size);
extern double summation(double* data, int size);
extern double* vector_fraction(double* data, int size, double denominator);
extern int comp_double (const void * elem1, const void * elem2);

extern double pvalue(double v, double* conddist, int size );

extern double fminsearch(pFMinSearch pfm, double* X0);

extern pFMinSearch fminsearch_new();
extern pFMinSearch fminsearch_new_with_eq(math_func eq, int Xsize, void* args);
extern void fminsearch_set_equation(pFMinSearch pfm,math_func eq, int Xsize, void* args);
extern void fminsearch_free(pFMinSearch pfm);
extern int fminsearch_min(pFMinSearch pfm, double* X0);
extern double* fminsearch_get_minX(pFMinSearch pfm);
extern double fminsearch_get_minF(pFMinSearch pfm);

typedef struct tagHistogram
{
	int nbins;
	int nsamples;
	double width;
	double max,min;

	double* point;
	unsigned int* count;
}Histogram;

typedef Histogram* pHistogram;

extern pHistogram histogram_new(double* data, int nsamples, int nbins);
extern void histogram_free(pHistogram phist);
extern void histogram_set_by_bin(pHistogram phist, double* data, int nsample, int nbins );
extern void histogram_set_by_unit(pHistogram phist, double* data, int nsample, double unit );
extern void histogram_set_sparse_data(pHistogram phist, double* data, int nsamples );
extern int histogram_get_count(pHistogram phist, double p);
extern double histogram_compare(pHistogram phist1, pHistogram phist2);
extern void histogram_print(pHistogram phist,FILE* fp);
extern int histogram_merge(pHistogram phist, pHistogram parg);
extern double histogram_check_fitness(pHistogram phist, double* args, double (*cdf)(double p, double* args) );
extern double histogram_get_prob(pHistogram phist, double p);
extern pHistogram histogram_load(char* file);
extern void histogram_set_with_preset_point(pHistogram phist, double* data, int nsamples, double* point, int nbins);

typedef struct tagANOVAResultElement
{
	int df;
	double* mean;	
	int* num;
	double vSS; 
	double MS;
	double F;
	double pvalue;
}ANOVAResultElement;

typedef ANOVAResultElement* pANOVAResultElement;

typedef struct tagANOVA
{
	int nways;
	int* ngrps;
	void* data;
	double gmean;
	pANOVAResultElement value;	
	ANOVAResultElement total;	
	ANOVAResultElement error;	
}ANOVA;

typedef ANOVA* pANOVA;

extern pANOVA anova_new(int nways, int* ngrps);
extern void anova(pANOVA panova, ... );
extern void anova1_run(pANOVA panova);
extern void anova2_run(pANOVA panova);
extern void anovan_run(pANOVA panova, va_list ap);
extern void anova_free(pANOVA panova);
extern void anova_print(pANOVA panova, char** name);
extern void anova_print_data(pANOVA panova);

extern double cmp_paired_t_test(double* grp1, double* grp2, int size);
extern double cmp_two_indep_chi2test(double* grp1, double* grp2, int size );
extern double cmp_two_indep_t_test(double* grp1, int size1, double* grp2, int size2 );

#ifdef	__cplusplus
}
#endif

#endif 
