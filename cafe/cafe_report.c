#include <mathfunc.h>
#include "cafe.h"
#include "viterbi.h"

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

void cafe_report_load_viterbi_pvalue(char* data, double** pvalues, int i, int nnodes)
{	
	char* next = &data[1];
	while( (next = index(next, ')' )) != NULL )
	{
		next++;
		if( next[0] == ',' )	
		{
			next[0] = '\t';
		}
		else if ( next[0] == ')' ) 
		{
			next[0] = '\0';
			break;
		}
	}
	double* array = cafe_report_load_data_double_pairs(&data[1],'\t');
	int j;
	for ( j = 0 ; j < nnodes - 1 ; j++ )
	{
		pvalues[j][i] = array[j];
	}
	memory_free(array);
	array = NULL;
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

int cafe_report_retrieve_data(const char* file, pCafeParam param)
{
	int i, j;
	family_size_range range;
	range.min = 0;
	range.max = 1;
	range.root_min = 1;
	range.root_max = 2;

	FILE* fp = fopen(file,"r");		
	if ( fp == NULL )
	{
		fprintf( stderr, "%s: Cannot open file: %s\n", __FUNCTION__, file );
		return -1;
	}
	pString pstr = string_new();
	int nnodes = 0;
	param->viterbi.expandRemainDecrease = (int**)memory_new(3,sizeof(int*));

	int bexist[2] = { 0 , 0 };
					
	while( file_read_line(pstr,fp) )
	{
		if ( pstr->buf[0] == '\'' ) 
		{
			pArrayList list = string_pchar_split(pstr->buf, '\t' );
			if ( list->size == 4 ) break;
			else if ( list->size == 6 )
			{
				bexist[0] = 1;	
				bexist[1] = 1;	
			}
			else
			{
				if ( (pstr->buf[1] & 0x4F) == 'C' )
				{
					bexist[0] = 1;	
				}
				else
				{
					bexist[1] = 1;	
				}
			}
			break;
		}
		char* data = index(pstr->buf, ':');
		if ( data == NULL )
		{
			fprintf( stderr, "Format error during loading cafe file: %s\n", pstr->buf );
			fclose(fp);
			return -2;
		}
		*data++ = '\0';			
		string_pchar_chomp(data);
		if ( strcasecmp( pstr->buf, "tree" ) == 0 )
		{
			param->pcafe = cafe_tree_new(data, &range, 0, 0);
			nnodes = param->pcafe->super.nlist->size;
		}
		else if ( strncasecmp( pstr->buf, "lambda tree", 10 ) == 0 )
		{
			param->lambda_tree = phylogeny_new( data,phylogeny_lambda_parse_func );
			tree_build_node_list(param->lambda_tree);
		}
		else if ( strncasecmp( pstr->buf, "lambda", 5 ) == 0  )
		{
			param->lambda = cafe_report_load_data_double_list(data, '\t', &param->num_lambdas );
		}
		else if ( strncasecmp( pstr->buf, "average", 7) == 0 )
		{
			param->viterbi.averageExpansion = cafe_report_load_data_double_pairs(data, '\t');
		}
		else if ( strncasecmp( pstr->buf, "expansion", 9 ) == 0 )
		{
			param->viterbi.expandRemainDecrease[0] = cafe_report_load_data_int_pairs(data,'\t');
		}
		else if ( strncasecmp( pstr->buf, "remain", 6 ) == 0 )
		{
			param->viterbi.expandRemainDecrease[1] = cafe_report_load_data_int_pairs(data,'\t');
		}
		else if ( strncasecmp( pstr->buf, "decrease", 8 ) == 0 )
		{
			param->viterbi.expandRemainDecrease[2] = cafe_report_load_data_int_pairs(data,'\t');
		}
	}
	param->param_set_func(param,param->lambda);
	pArrayList plines = arraylist_new(11000);
	int num_families;
	for ( num_families  = 0 ;  file_read_line( pstr, fp ); num_families++ )
	{
		char* line = (char*)memory_new(pstr->length+1,sizeof(char));
		strcpy( line, pstr->buf );
		arraylist_add(plines, line); 				
	}

	viterbi_parameters_init(&param->viterbi, nnodes, num_families);

	param->pfamily = (pCafeFamily)memory_new(1, sizeof(CafeFamily) );
	pCafeFamily pcf = param->pfamily;
	pcf->flist = arraylist_new(11000);
	pcf->num_species = (nnodes + 1)/2;
	pcf->species = (char**)memory_new( pcf->num_species, sizeof(char*));
	pcf->index = (int*)memory_new( pcf->num_species ,sizeof(int) );
	pArrayList nlist = param->pcafe->super.nlist;	

	
	for ( i = j = 0 ; i < nnodes ; i+=2, j++ )
	{
		char* name = ((pPhylogenyNode)nlist->array[i])->name;
		pcf->index[j] =  i;
		pcf->species[j] = (char*)memory_new( strlen(name) + 1, sizeof(char));
		strcpy( pcf->species[j], name );
	}	

	if ( bexist[0] )
	{
		param->viterbi.cutPvalues = (double**)memory_new_2dim(nnodes, num_families, sizeof(double) );
	}	
	if ( bexist[1] )
	{
		param->likelihoodRatios = (double**)memory_new_2dim(nnodes, num_families, sizeof(double) );		
	}
	int max_size = 0;
	for ( i = 0 ; i < plines->size ; i++ )
	{
		pArrayList data = string_pchar_split( (char*)plines->array[i], '\t' );
		pCafeFamilyItem pitem = (pCafeFamilyItem)memory_new(1,sizeof(CafeFamilyItem));
		pitem->desc = NULL;
		pitem->id = (char*)data->array[0];		
		pitem->count = (int*)memory_new(pcf->num_species, sizeof(int));
		pitem->maxlh = -1;
		pitem->ref = -1;
		pitem->lambda = NULL;
		arraylist_add(pcf->flist,pitem);
		sscanf((char*)data->array[2], "%lf", &param->viterbi.maximumPvalues[i]);

		pCafeTree ptree = cafe_tree_new( (char*)data->array[1], &range, 0, 0 );
		pArrayList nlist = ptree->super.nlist;

		for ( j = 0 ; j < nlist->size ; j+=2 )
		{
			pCafeNode pcnode = (pCafeNode)nlist->array[j];
			pitem->count[j/2] = pcnode->familysize;
			if ( max_size < pcnode->familysize )
			{
				max_size = pcnode->familysize;
			}
		}
		for ( j = 1 ; j < nlist->size ; j+= 2 )
		{
			pCafeNode pcnode = (pCafeNode)nlist->array[j];
			param->viterbi.viterbiNodeFamilysizes[j/2][i] = pcnode->familysize;
		}

		cafe_report_load_viterbi_pvalue((char*)data->array[3], param->viterbi.viterbiPvalues, i, nnodes );
		if ( bexist[0] )
		{
			cafe_report_load_bc_or_lhr_list((char*)data->array[4], 
					    param->viterbi.cutPvalues, i, nnodes );
			param->viterbi.cutPvalues[param->pcafe->super.root->id][i] = -1.0;
		}
		if ( bexist[1] )
		{
			int idx =  bexist[0] ? 5 : 4;
			cafe_report_load_bc_or_lhr_list((char*)data->array[idx], 
					    param->likelihoodRatios, i, nnodes );
		}
		cafe_tree_free(ptree);
	}
	
	init_family_size(&range, max_size);

	cafe_tree_set_parameters(param->pcafe, &range, param->lambda[0] );
	arraylist_free(plines, free);
	return 0;
}
