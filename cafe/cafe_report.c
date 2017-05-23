#include <mathfunc.h>
#include "cafe.h"

/**************************************************************************
 * Retrieve
**************************************************************************/

extern void phylogeny_lambda_parse_func(pTree ptree, pTreeNode ptnode);

double* cafe_report_load_data_double_list(char* data, int delimiter, int* num)
{
	pArrayList list = string_pchar_split(data, delimiter);	
	double* array = (double*)memory_new(list->size,sizeof(double));
	int i;
	if( num ) *num = list->size;
	for ( i = 0 ; i < list->size ; i++ )
	{
		array[i] = -1;
		sscanf( (char*)list->array[i], "%lf", &array[i] );
	}
	arraylist_free(list, free);
	return array;
}

double* cafe_report_load_data_double_pairs(char* data, int delimiter)
{
	pArrayList list = string_pchar_split(data, delimiter);	
	double* array = (double*)memory_new(list->size*2,sizeof(double));
	int i;
	for ( i = 0 ; i < list->size ; i++ )
	{
		array[2*i] = -1;
		array[2*i+1] = -1;
		sscanf( (char*)list->array[i], "(%lf,%lf)", &array[i*2], &array[i*2+1] );
	}
	arraylist_free(list, free);
	return array;
}

int* cafe_report_load_data_int_pairs(char* data, int delimiter)
{
	pArrayList list = string_pchar_split(data, delimiter);	
	int* array = (int*)memory_new(list->size*2,sizeof(int));
	int i;
	for ( i = 0 ; i < list->size ; i++ )
	{
		array[i*2] = -1;
		array[i*2+1] = -1;
		sscanf( (char*)list->array[i], "(%d,%d)", &array[i*2], &array[i*2+1] );
	}
	arraylist_free(list, free);
	return array;
}

void cafe_report_load_bc_or_lhr_list(char* data, double** pvalues, int i, int nnodes)
{
	size_t len = strlen(data);
	data[len-1] = '\0';
	double* array = cafe_report_load_data_double_list(&data[1],',', NULL);
	int j;
	for ( j = 0 ; j < nnodes; j++ )
	{
		pvalues[j][i] = array[j];
	}
	memory_free(array);
	array = NULL;
}

size_t file_read_line(pString pstr, FILE* fp)
{
	char buf[STRING_BUF_SIZE];
	string_reset(pstr);
	while (fgets(buf, STRING_BUF_SIZE, fp))
	{
		string_add(pstr, buf);
		size_t len = strlen(buf);
		if (buf[len - 1] == '\n') break;
	}
	return pstr->length;
}

