#include<mathfunc.h>
#include<memalloc.h>


/***********************************************************************
 * Histogram
 ***********************************************************************/

pHistogram histogram_new(double* data, int nsamples, int nbins)
{
	pHistogram phist = (pHistogram) memory_new(1, sizeof(Histogram));
	if ( data == NULL ) 
	{
		return phist;
	}
	phist->nbins = nbins;
	if ( nbins > 0 )
	{
		histogram_set_by_bin(phist,data,nsamples,nbins);
	}
	else if ( nbins == -1 )
	{
		histogram_set_sparse_data(phist,data,nsamples);
	}
	return phist;
}

pHistogram histogram_load(char* file)
{
	FILE* fp;
	if ( (fp = fopen(file,"r"))  == NULL )
	{
		perror(file);	
		return NULL;
	}
	pHistogram phist = (pHistogram) memory_new(1, sizeof(Histogram));
	fscanf(fp,"MIN: %lf ~ MAX: %lf\n", &phist->min , &phist->max );	
	fscanf(fp,"BIN: %d\n", &phist->nbins );
	fscanf(fp,"# Samples: %d\n", &phist->nsamples);
	phist->point = (double*) memory_new(phist->nbins, sizeof(double) );
	phist->count = (unsigned int*) memory_new(phist->nbins, sizeof(unsigned int));
	int i;
	float f;
	for ( i = 0 ; i < phist->nbins ; i++ )
	{
		fscanf(fp,"%lf\t%d\t%g\n", &phist->point[i],&phist->count[i],&f);
	}
	fclose(fp);
	return phist;
}

void histogram_reset(pHistogram phist)
{
	if (phist->point) memory_free(phist->point);
	phist->point = NULL;
	if (phist->count) memory_free(phist->count);
	phist->count = NULL;
}

void histogram_free(pHistogram phist)
{
	histogram_reset(phist);
	memory_free( phist );
	phist = NULL;
}


void __max_and_min(double* data, int size, double* max, double* min)
{
	int i;
	*max = *min = data[0];
	for ( i = 1 ; i < size ; i++ )
	{
		if ( data[i] > *max ) 
		{
			*max = data[i];
		}
		else if ( data[i] < *min )
		{
			*min = data[i];
		}
	}
}

int __histogram_index(pHistogram phist, double p)
{
	if ( p == phist->max ) return phist->nbins-1;
	else if ( p == phist->min ) return 0;

	if ( phist->point == NULL )
	{
		return (int)((p - phist->min)/phist->width);
	}
	int i;
	if ( phist->width == 0 )
	{
		if ( phist->max < p || phist->min > p ) return -1;
		int bottom = 1;
		int top = phist->nbins - 2;
		int cur = phist->nbins/2;
		while(1)
		{
			if ( phist->point[cur] <= p && phist->point[cur+1] > p ) return cur;
			else if ( phist->point[cur] > p && phist->point[cur-1] <= p ) return cur-1;
			if ( phist->point[cur] > p )
			{
				top = cur-1; 
				cur = bottom + (cur - bottom)/2;
			}
			else
			{
				bottom = cur + 1;
				cur = bottom + (top-cur)/2;
			}
		}
	}
	else
	{
		for ( i = 0 ; i < phist->nbins ; i++ )
		{
			if ( phist->point[i] == p  ) return i;
		}
	}
	return -1;
}

void histogram_set_by_bin(pHistogram phist, double* data, int nsamples, int nbins )
{
	int i;
	if ( phist->point ) 
	{
		memory_free(phist->point);
		phist->point = NULL;
	}
	if ( phist->nbins != nbins || phist->count == NULL )
	{
		if ( phist->count ) memory_free(phist->count);
		phist->nbins = nbins;
		phist->count = (unsigned int*) memory_new( nbins, sizeof(unsigned int) );
	}
	else
	{
		memset( phist->count , 0, sizeof(unsigned int)*nbins );
	}
	__max_and_min(data, nsamples, &phist->max, &phist->min );
	phist->nsamples = nsamples;
	phist->width = (phist->max - phist->min)/nbins;
	for ( i = 0 ; i < nsamples; i++ )
	{
		phist->count[__histogram_index(phist,data[i])]++;
	}
}

void histogram_set_by_unit(pHistogram phist, double* data, int nsamples, double unit )
{
	int i;
	__max_and_min(data, nsamples, &phist->max, &phist->min );
	phist->nsamples = nsamples;
	phist->width = unit;
	phist->max = ceil(phist->max/unit)*unit;
	phist->min = floor(phist->min/unit)*unit;
	int nbins = (int)((phist->max - phist->min+unit)/unit);

	if( phist->point ) 
	{
		memory_free(phist->point);
		phist->point = NULL;
	}
	if ( phist->nbins != nbins || phist->count == NULL )
	{
		if ( phist->count ) memory_free(phist->count);
		phist->nbins = nbins;
		phist->count = (unsigned int*) memory_new( phist->nbins, sizeof(unsigned int) );
	}
	else
	{
		memset( phist->count , 0, sizeof(unsigned int)*nbins );
	}

	for ( i = 0 ; i < nsamples ; i++ )
	{
		phist->count[__histogram_index(phist,data[i])]++;
	}
}

void __qsort_double_with_index(double* list, int* idx, int left, int right);

int* qsort_double_with_index(double* list, int size)
{
	int i;
	int* index = (int*) memory_new(size,sizeof(int));
	for( i = 0 ; i < size; i++ ) index[i] = i;
	__qsort_double_with_index(list, index, 0, size-1);
	return index;	
}

void histogram_set_with_preset_point(pHistogram phist, double* data, int nsamples, double* point, int nbins)
{
	phist->width = 0;
	phist->nbins = nbins;
	phist->min = point[0];
	phist->max = point[nbins-1];
	phist->point = (double*) memory_new(nbins,sizeof(double));
	memcpy( phist->point, point, sizeof(double)*nbins);
	phist->count = (unsigned int*) memory_new(nbins,sizeof(unsigned int));
	phist->nsamples = nsamples;
	int i;
	for( i = 0 ; i < nsamples; i++ )
	{
		int idx = __histogram_index(phist,data[i]);
		if ( idx == -1 ) continue;
		phist->count[idx]++;
	}
}

void histogram_set_sparse_data(pHistogram phist, double* data, int nsamples )
{
	int size = 100;
	int i,j, idx;
	histogram_reset(phist);
	phist->nsamples = nsamples;
	phist->nbins = 0;
	phist->width = -1;
	phist->point = (double*) memory_new(size,sizeof(double));	
	phist->count = (unsigned int*) memory_new(size,sizeof(unsigned int));
	for ( i = 0 ; i < nsamples ; i++ )
	{
		idx = -1;
		for ( j = 0 ; j < phist->nbins ; j++ )
		{
			if ( phist->point[j] == data[i] ) 
			{
				idx = j;
				break;
			}
		}
		if ( idx != -1 )
		{
			phist->count[idx]++;
		}
		else
		{
			if ( size == 0 )
			{
				phist->point = (double*)memory_realloc(phist->point,phist->nbins+100,sizeof(double));	
				phist->count = (unsigned int*)memory_realloc(phist->count,phist->nbins+100,sizeof(unsigned int));
				size = 100;
			}
			size--;
			phist->point[phist->nbins] = data[i];
			phist->count[phist->nbins]++;
			phist->nbins++;
		}
	}
	unsigned int* index = (unsigned int*)qsort_double_with_index(phist->point,phist->nbins);
	for( i = 0 ; i < phist->nbins ; i++ )
	{
		index[i] = phist->count[index[i]];
	}
	memory_free(phist->count);
	phist->count = NULL;
	phist->count = index;
	phist->max = phist->point[phist->nbins-1];
	phist->min = phist->point[0];
}

int histogram_get_count(pHistogram phist, double p)
{
	if ( phist->max < p || phist->min > p ) return 0;
	int idx = __histogram_index(phist,p);
	return idx == -1 ? 0 : phist->count[idx];
}

double histogram_get_prob(pHistogram phist, double p)
{
	return histogram_get_count(phist,p)/(double)phist->nsamples;
}

int histogram_merge(pHistogram phist, pHistogram parg)
{
	if ( phist->nsamples == 0 ) 
	{
		memcpy( phist, parg, sizeof(Histogram));
		phist->count = (unsigned int*) memory_new( phist->nbins,sizeof(unsigned int));
		memcpy( phist->count, parg->count, sizeof(unsigned int)*phist->nbins );
		if ( parg->point )
		{
			phist->point = (double*) memory_new(phist->nbins,sizeof(double));
			memcpy( phist->point, parg->point, sizeof(double)*phist->nbins );
		}
		return 0;
	}
	if ( parg->point )
	{
		if ( phist->point == NULL )
		{
			fprintf(stderr, "ERROR(histogram_merge): histogram must have the same type of histogram: sparse or regular\n");
			return -1;
		}
		int i, idx;	
		phist->nsamples += parg->nsamples;
		int addsize = 0;
		for ( i = 0 ; i < parg->nbins ; i++ )
		{
			idx = __histogram_index(phist, parg->point[i]);
			if ( idx == -1 )
			{
				addsize++;	
			}
			else
			{
				phist->count[idx] += parg->count[i];
			}
		}
		if ( addsize ==  0) return 0;
		int k = phist->nbins;
		phist->nbins += addsize;
		phist->point = (double*)memory_realloc(phist->point,phist->nbins,sizeof(double));	
		phist->count = (unsigned int*)memory_realloc(phist->count,phist->nbins,sizeof(unsigned int));
		for ( i = 0 ; i < parg->nbins ; i++ )
		{
			idx = __histogram_index(phist, parg->point[i]);
			if ( idx == -1 || idx >= k)
			{
				phist->point[k] = parg->point[i];
				phist->count[k] = parg->count[i];
				k++;
			}
		}
		unsigned int* index = (unsigned int*)qsort_double_with_index(phist->point,phist->nbins);
		for( i = 0 ; i < phist->nbins ; i++ )
		{
			index[i] = phist->count[index[i]];
		}
		memory_free(phist->count);
		phist->count = NULL;
		phist->count = index;
		phist->max = phist->point[phist->nbins-1];
		phist->min = phist->point[0];
	}
	else
	{

	}
	return 0;
}

void histogram_print(pHistogram phist, FILE* fp)
{
	int i;
	if( fp == NULL ) fp = stderr;
	fprintf(fp,"MIN: %f ~ MAX: %f\n", phist->min, phist->max  );
	fprintf(fp,"BIN: %d\n", phist->nbins );
	fprintf(fp,"# Samples: %d\n", phist->nsamples );
	if ( phist->point )
	{
		for ( i = 0 ; i < phist->nbins; i++ )
		{
			fprintf(fp,"%g\t%d\t%g\n", phist->point[i], phist->count[i], ((double)phist->count[i])/phist->nsamples );
		}
	}
	else
	{
		for ( i = 0 ; i < phist->nbins ; i++ )
		{
			fprintf(fp,"%g\t%d\t%g\n", phist->min+phist->width*i, phist->count[i],
								((double)phist->count[i])/phist->nsamples );
		}
	}
}

double histogram_compare(pHistogram phist1, pHistogram phist2)
{
	return 0;
}

double histogram_check_fitness(pHistogram phist, double* args, double (*cdf)(double p, double* args) )
{
	if ( phist->nbins < 5 )
	{
	//	fprintf(stderr, "The number of bins must be greater than or equal to 5\n" );
		return -1;
	}
	int i;
	double wing[2] = { 0 , 0 };
	double* e = (double*) memory_new(phist->nbins, sizeof(double));
	double chi2 = 0;
	for ( i = 0 ; i < phist->nbins ; i++ )
	{
		e[i] = cdf( i*phist->width + phist->min, args);
	}

	wing[0] = e[0];
	wing[1] = cdf(phist->max, args);
	for ( i = 0 ; i < phist->nbins - 1 ; i++ )
	{
		e[i] = e[i+1] - e[i];
	}
	e[phist->nbins-1] = wing[1] - e[phist->nbins-1];
	wing[1] = 1 - wing[1];

	/*
	fprintf(stderr, "%f %f\n", wing[0], wing[0] );
	double sum  = wing[0];
	for ( i = 0 ; i < phist->nbins ; i++ )
	{
		sum += e[i];
		fprintf(stderr, "%f, %f\n", e[i], sum );
	}
	fprintf(stderr, "%f %f\n", wing[1], wing[1]+sum );
	*/

	for ( i = 0 ; i < 2 ; i++ )
	{
		chi2 += wing[i]*phist->nsamples;
	}
	for ( i = 0 ; i < phist->nbins ; i++ )
	{
		e[i] *= phist->nsamples;
		chi2 += ( phist->count[i] - e[i] ) * ( phist->count[i] - e[i] ) / e[i];
	}
//	fprintf(stderr, "%f, %f, %d\n", chi2, chi2cdf(chi2,phist->nbins-1), phist->nbins );
	memory_free(e);
	e = NULL;
	return 1-chi2cdf(chi2,phist->nbins-1-2);
}
