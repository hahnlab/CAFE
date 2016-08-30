#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<utils.h>
#include<mathfunc.h>
#include<memalloc.h>


double cmp_paired_t_test(double* grp1, double* grp2, int size)
{
	int i;
	double mean = 0;
	double std = 0;
	for ( i = 0 ; i < size; i++ )
	{
		mean += grp1[i] - grp2[i];
	}
	mean /= size;
	for ( i = 0 ; i < size ; i++ )
	{
		std += pow( grp1[i] - grp2[i] - mean,2);
	}
	std = sqrt(std/(size*(size-1)));
	fprintf(stderr, "mean: %f, std: %f\n", mean , std );
	return mean > 0 ? 1 - tcdf(mean/std, size-1) : tcdf(mean/std,size-1);	
}

double cmp_two_indep_t_test(double* grp1, int size1, double* grp2, int size2 )
{
	int i, j;
	double* table[2] = { grp1, grp2 };
	int n[2] = { size1, size2 };
	double y[2], v[2];
	for( i = 0 ; i < 2 ; i++ )
	{	
		y[i] = 0;
		for( j = 0 ; j < n[i] ; j++ )
		{
			y[i] += table[i][j];
		}
		y[i] /= n[i];
	}
	for( i = 0 ; i < 2 ; i++ )
	{	
		v[i] = 0;
		for( j = 0 ; j < n[i] ; j++ )
		{
			v[i] += pow(table[i][j] - y[i],2);
		}
		v[i] /= (n[i]-1);
	}
	int df = n[0] + n[1] - 2;
	double s = sqrt((( n[0] - 1 ) * v[0]+ (n[1]- 1 ) * v[1])/df);
	s *= sqrt(1.0/n[0]+ 1.0/n[1]);
	return 1 - tcdf((y[0]-y[1])/s,df);
}


/******************************************************************************
 * Compare two independent sample  : H0 u11 = u21, u21 = u22, u31= u32 ...
 * ***************************************************************************/
double cmp_two_indep_chi2test(double* grp1, double* grp2, int size )
{
	double p[2] = { summation(grp1,size), summation(grp2,size) };	
	double* table[2] = { grp1, grp2 };
	double total = p[0] + p[1];
	p[0] /= total; p[1] /= total;
	int i, j;
	double chi2 = 0;
	for ( j = 0 ; j < size ; j++ )
	{
		double csum = table[0][j] + table[1][j];
		for ( i = 0 ; i < 2 ; i++ )
		{
			double e = csum * p[i];
			chi2 += ( table[i][j] - e ) * ( table[i][j] - e ) / e;
		}
	}
	return 1 - chi2cdf( chi2, size - 1 );
}	

/******************************************************************************
 * ANOVA
 * ***************************************************************************/

pANOVA anova_new(int nways, int* ngrps)
{
	pANOVA panova = (pANOVA) memory_new(1,sizeof(ANOVA));	

	panova->nways = nways;
	nways = abs(nways);
	panova->ngrps = (int*) memory_new(nways,sizeof(int));
	panova->gmean = 0;

	int i;
	int nvalues = (1<<abs(nways)) - 1;
	double*** data;
	
	switch( panova->nways )
	{
		case 1:
			panova->data = (double**) memory_new( ngrps[0], sizeof(double*) );
			panova->ngrps[0] = ngrps[0];
		   	break;
		case 2:
			panova->data = (double***) memory_new( ngrps[0], sizeof(double**) );
			data = (double***)panova->data;
			for( i = 0 ; i < ngrps[0] ; i++ )
			{
				data[i] = (double**) memory_new( ngrps[1], sizeof(double*));
			}
			memcpy( panova->ngrps, ngrps, sizeof(int)*nways );
		   	break;
		default:
			panova->data = NULL;
			memcpy( panova->ngrps, ngrps, sizeof(int)*nways );
	}
	panova->value = (pANOVAResultElement) memory_new( nvalues, sizeof(ANOVAResultElement) );
	return panova;
}

void anova_free(pANOVA panova)
{
	int i;
	int size = (1<<panova->nways)-1;
	for ( i = 0 ; i < size ; i++ )
	{
		if ( panova->value[i].mean ) memory_free(panova->value[i].mean);
		panova->value[i].mean = NULL;
		if ( panova->value[i].num ) memory_free(panova->value[i].num );
		panova->value[i].num = NULL;
	}
	double*** data;
	switch( panova->nways )
	{
		case 2:
			data = (double***)panova->data;
			for ( i = 0 ; i < panova->ngrps[0] ; i++ )
			{
				memory_free(data[i]);
				data[i] = NULL;
			}
		case 1:
			memory_free(panova->data);
			panova->data = NULL;
			break;
	}
	memory_free(panova->ngrps);
	panova->ngrps = NULL;
	memory_free(panova->value);
	panova->value = NULL;
	memory_free(panova);
	panova = NULL;
}

void anova_degree_of_freedom(pANOVA panova, va_list ap1)
{
	int i, j;
	double** data1;
	double*** data2;
  va_list ap;
  va_copy(ap, ap1);
	switch( panova->nways )
	{
		case 1:
			data1 = (double**)panova->data;
			panova->total.df = 0;
			for ( i = 0 ; i < panova->ngrps[0] ; i++ )
			{
				panova->total.df +=  (int)data1[i][0];
			}	
			panova->total.df--;
			panova->value[0].df = panova->ngrps[0] - 1;
			panova->error.df = panova->total.df - panova->value[0].df; 
			break;
		case 2:
			data2 = (double***)panova->data;
			panova->total.df = 0;
			for ( i = 0 ; i < panova->ngrps[0]	; i++ )
			{
				for ( j = 0 ; j < panova->ngrps[1] ; j++ )
				{
					panova->total.df += (int)data2[i][j][0];
				}
			}
			panova->total.df--;
			panova->value[0].df = panova->ngrps[0] - 1;
			panova->value[1].df = panova->ngrps[1] - 1;
			panova->value[2].df = panova->value[0].df * panova->value[1].df;
			panova->error.df = panova->total.df - ( panova->value[0].df 
							   + panova->value[1].df 
				               + panova->value[2].df); 
			break;
		default:
			break;
	}
  va_end(ap);
}

void anova_group_error_and_stat(pANOVA panova)
{
	int i, j;
	int nways = abs(panova->nways);
	int nvalues = (1<<nways)-1;

	for ( i = 0 ; i < nways ; i++ )
	{
		for ( j = 0 ; j < panova->ngrps[i]; j++ )
		{
			panova->value[i].vSS += panova->value[i].num[j] * pow( panova->value[i].mean[j] - panova->gmean, 2);
		}
	}

	panova->error.vSS = panova->total.vSS;
	for ( i = 0 ; i < nvalues ; i++ )
	{
		panova->error.vSS -= panova->value[i].vSS;
		panova->value[i].MS = panova->value[i].vSS/panova->value[i].df;
	}
	panova->error.MS = panova->error.vSS/panova->error.df;

	for ( i = 0 ; i < nvalues ; i++ )
	{
		panova->value[i].F = panova->value[i].MS/panova->error.MS;
		panova->value[i].pvalue = 1 - fcdf( panova->value[i].F, 
			                            panova->value[i].df,
			                            panova->error.df );
	}
}

void anova(pANOVA panova, ... )
{
	va_list ap;
	va_start(ap, panova);
	anova_degree_of_freedom(panova,ap);
	switch( panova->nways )
	{
		case 1: anova1_run(panova); break;
		case 2: anova2_run(panova); break;
		default: break;
	}
	anova_group_error_and_stat(panova);
	va_end(ap);
}

void anova1_run(pANOVA panova)
{
	int i, k;
	double** table = (double**)panova->data;

// Mean
	int cnt = 0;
	panova->value[0].num = (int*) memory_new(panova->ngrps[0], sizeof(int));
	panova->value[0].mean = (double*) memory_new(panova->ngrps[0], sizeof(double));
	panova->gmean = 0;
	for ( i = 0 ; i < panova->ngrps[0] ; i++ )
	{
		for ( k = 1 ; k <= table[i][0] ; k++ )
		{
			panova->value[0].mean[i] += table[i][k];
			panova->gmean += table[i][k];
		}
		panova->value[0].mean[i] /= table[i][0];
		panova->value[0].num[i] = (int)table[i][0];
		cnt += (int)table[i][0];
	}
	panova->gmean /= cnt;

// vSS
	for ( i = 0 ; i < panova->ngrps[0]; i++ )
	{
		for ( k = 1 ; k <= table[i][0] ; k++ )
		{
			panova->total.vSS += pow(table[i][k] - panova->gmean ,2);
		}		
	}
}

void anova2_run(pANOVA panova)
{
	int i,j, k;
	double*** table = (double***)panova->data;

// Means
	for ( i = 0 ; i < 2 ; i++ )
	{
		panova->value[i].mean = (double*) memory_new(panova->ngrps[i],sizeof(double));
		panova->value[i].num = (int*) memory_new(panova->ngrps[i],sizeof(int));
	}
	panova->gmean = 0;
	panova->value[2].mean = (double*) memory_new(panova->ngrps[0]*panova->ngrps[1],sizeof(double));
	
	int cnt = 0;
	for ( i = 0 ; i < panova->ngrps[0] ; i++ )
	{
		int base = i * panova->ngrps[1];
		int subcnt = 0;
		for ( j = 0 ; j < panova->ngrps[1] ; j++ )
		{
			for ( k = 1 ; k <= table[i][j][0] ; k++ )
			{
				panova->value[2].mean[base+j] += table[i][j][k];
				panova->value[0].mean[i] += table[i][j][k];
				panova->value[1].mean[j] += table[i][j][k];
				panova->gmean += table[i][j][k];
			}
			panova->value[2].mean[base+j] /= table[i][j][0];
			cnt += (int)table[i][j][0];
			subcnt += (int)table[i][j][0];
		}
		panova->value[0].num[i] = subcnt;
		panova->value[0].mean[i] /= subcnt;
	}
	panova->gmean /= cnt;
	for ( j = 0 ; j < panova->ngrps[1] ; j++ )
	{
		int cnt= 0;
		for ( i = 0 ; i < panova->ngrps[0] ; i++ )
		{
			cnt += (int)table[i][j][0];
		}
		panova->value[1].mean[j] /= cnt;
		panova->value[1].num[j] = cnt;
	}

// vSS
	for ( i = 0 ; i < panova->ngrps[0] ; i++ )
	{
		int base = i * panova->ngrps[1];
		for ( j = 0 ; j < panova->ngrps[1] ; j++ )
		{
			panova->value[2].vSS += table[i][j][0] 
				    * pow( panova->value[2].mean[base+j] - panova->value[0].mean[i] 
						- panova->value[1].mean[j] + panova->gmean, 2);
			for ( k = 1 ; k <= table[i][j][0] ; k++ )
			{
				panova->total.vSS += pow(table[i][j][k] - panova->gmean,2);
			}
		}
	}
}

void anova_print(pANOVA panova, char** name)
{
	int i;
	int nvalues = (1<<abs(panova->nways))-1;
	fprintf(stderr, "Source\tSS\tdf\tMS\tF\tp-value\n");
	for ( i = 0 ; i < nvalues ; i++ )
	{
		if ( name ) printf("%s", name[i] );
		else printf("%d", i+1);
		fprintf(stderr, "\t%f\t%d\t%f\t%f\t%f\n", panova->value[i].vSS,
									panova->value[i].df, panova->value[i].MS,
									panova->value[i].F, panova->value[i].pvalue );
	}
	fprintf(stderr, "Error\t%f\t%d\t%f\n", panova->error.vSS, panova->error.df, panova->error.MS );
	fprintf(stderr, "Total\t%f\t%d\n", panova->total.vSS, panova->total.df);
}

void anova_print_data(pANOVA panova)
{
	if ( panova->nways == 1 )
	{
		int i, j;
		double** data = (double**)panova->data;
		double max = data[0][0];
		for ( i = 1 ; i < panova->ngrps[0] ; i++ )
		{
			max = data[i][0] > max ? data[i][0] : max;
		}
		for ( j = 1 ; j <= max ; j++ )
		{
			fprintf(stderr, "%d", j );
			for ( i = 0 ; i < panova->ngrps[0] ; i++ )
			{
				data[i][0] >= j ?  printf("\t%f", data[i][j]) : printf("\t");
			}
			fprintf(stderr, "\n");
		}
	}
}
