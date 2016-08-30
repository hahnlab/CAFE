#include<cafe.h>
#include<io.h>

/**************************************************************************
 * Save
**************************************************************************/

void cafe_report_set_viterbi(pCafeParam param, int i)
{
	int j;
	cafe_family_set_size(param->pfamily,i, param->pcafe);
	pArrayList nlist = param->pcafe->super.nlist;
	for ( j = 1 ; j < nlist->size ; j+=2 )
	{
		pCafeNode pcnode = (pCafeNode)nlist->array[j];
		pcnode->familysize = param->viterbiNodeFamilysizes[j/2][i];
	}
}

void lambda_tree_string(pString pstr, pPhylogenyNode pnode )
{
	if ( pnode->taxaid != -1 ) 
	{
		string_fadd( pstr, "%d", pnode->taxaid+1 );
	}
}

void cafe_report_text(pCafeParam param)
{
	int i, b;
	pArrayList pflist = param->pfamily->flist;
	pArrayList nlist = param->pcafe->super.nlist;

	fprintf(param->fout, "Tree:");
	pString pstr = phylogeny_string((pTree)param->pcafe, NULL);	
	fprintf(param->fout,"%s\n", pstr->buf );
	string_free(pstr);

	fprintf(param->fout, "Lambda:");
	for ( i = 0 ;  i < param->num_lambdas ; i++ )
	{
		fprintf(param->fout, "\t%g", param->lambda[i] );
	}
	fprintf( param->fout, "\n");
	if ( param->lambda_tree )
	{
		pString pstr = phylogeny_string_newick(param->lambda_tree, lambda_tree_string, 0 );	
		fprintf(param->fout, "Lambda tree:\t%s\n", pstr->buf );
		string_free(pstr);
	}

	fprintf(param->fout, "# IDs of nodes:");
//	pstr = cafe_tree_string_with_id(param->pcafe);
	void cafe_tree_string_id(pString pstr, pPhylogenyNode pnode);
	pstr = phylogeny_string_newick((pTree)param->pcafe, cafe_tree_string_id ,PS_SKIP_BL );
	fprintf(param->fout,"%s\n", pstr->buf );
	string_free(pstr);

	int nnodes = nlist->size/2;
    fprintf(param->fout,"# Output format for: ' Average Expansion', 'Expansions', 'No Change', 'Contractions', and 'Branch-specific P-values' = (node ID, node ID): ");
  	for ( b = 1 ; b < nlist->size ; b+=2 )
	{
		pTreeNode pnode = (pTreeNode)nlist->array[b];
		pTreeNode child[2] = { tree_get_child(pnode,0), tree_get_child(pnode,1) };
		fprintf(param->fout,"(%d,%d) ", child[0]->id, child[1]->id );
	}
    fprintf(param->fout,"\n");
	fprintf(param->fout,"# Output format for 'Branch cutting P-values' and 'Likelihood Ratio Test': (%d", 0);
	for ( i = 1 ; i < nlist->size; i++ )
	{
		fprintf(param->fout, ", %d", i );
	}
	fprintf(param->fout, ")\n");
/*
	fprintf(param->fout,"# Branch located between parent and child for expansion/remain/decrease and Node P-value\n");
	fprintf(param->fout,"# Format: parant id:(left child id,right child id)\n");
	fprintf(param->fout,"# ");
	for ( b = 1 ; b < nlist->size ; b+=2 )
	{
		pTreeNode pnode = (pTreeNode)nlist->array[b];
		pTreeNode child[2] = { tree_get_child(pnode,0), tree_get_child(pnode,1) };
		fprintf(param->fout,"%d:(%d,%d) ", b, child[0]->id, child[1]->id );
	}
	fprintf(param->fout,"\n");
*/
	fprintf(param->fout, "Average Expansion:");
	for ( b = 0 ; b < nnodes ; b++ )
	{
		fprintf(param->fout,"\t(%lf,%lf)", param->averageExpansion[2*b], 
						param->averageExpansion[2*b+1] );
	}
	fprintf(param->fout,"\nExpansion :");
	for ( b = 0 ; b < nnodes ; b++ )
	{
		fprintf(param->fout,"\t(%d,%d)", param->expandRemainDecrease[0][2*b], 
						param->expandRemainDecrease[0][2*b+1] );
	}
	fprintf(param->fout,"\nRemain :");
	for ( b = 0 ; b < nnodes ; b++ )
	{
		fprintf(param->fout,"\t(%d,%d)", param->expandRemainDecrease[1][2*b], 
						param->expandRemainDecrease[1][2*b+1] );
	}
	fprintf(param->fout,"\nDecrease:");
	for ( b = 0 ; b < nnodes ; b++ )
	{
		fprintf(param->fout,"\t(%d,%d)", param->expandRemainDecrease[2][2*b], 
						param->expandRemainDecrease[2][2*b+1] );
	}
	fprintf(param->fout,"\n");
/*
	fprintf(param->fout,"# cut P-value and Likelihodd ratio : ");
	for ( i = 0 ; i < nlist->size ; i++ )
	{
		pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
		if ( pnode->name ) 
		{
			fprintf(param->fout,"   %s:", pnode->name );
		}
		else if ( tree_is_root( (pTree)param->pcafe, (pTreeNode)pnode ) )
		{
			fprintf(param->fout,"    root:" );
		}
		else
		{
			fprintf(param->fout,"   ");
		}
		fprintf(param->fout,"%d", pnode->super.id );
	}
	fprintf(param->fout,"\n");
*/
	fprintf(param->fout, "'ID'\t'Newick'\t'Family-wide P-value'\t'Viterbi P-values'");
	if ( param->cutPvalues )
	{
		fprintf(param->fout, "\t'cut P-value'");
	}
	if ( param->likelihoodRatios )
	{
		fprintf(param->fout, "\t'Likelihood Ratio'");
	}
	fprintf(param->fout, "\n");
	for ( i = 0 ; i < pflist->size ; i++ )
	{
		fprintf(param->fout,"%s\t", ((pCafeFamilyItem)pflist->array[i])->id );
		cafe_report_set_viterbi(param,i);
		pString pstr = cafe_tree_string(param->pcafe);	
		fprintf(param->fout,"%s\t", pstr->buf );
		string_free(pstr);
		fprintf(param->fout,"%lf\t(", param->maximumPvalues[i] );
		for ( b = 0 ; b < nnodes ; b++ )
		{
			if ( param->viterbiPvalues[2*b][i] == -1 )
			{
				fprintf(param->fout,"(-,-)");
			}
			else
			{
				fprintf(param->fout,"(%lf,%lf)", param->viterbiPvalues[2*b][i], param->viterbiPvalues[2*b+1][i]);
			}
			if ( b < nnodes - 1 ) fprintf( param->fout, "," );
		}
		fprintf(param->fout,")\t" );

		if ( param->cutPvalues )
		{
			fprintf(param->fout,"(" );
			for ( b = 0 ; b < nlist->size; b++ )
			{
				if ( param->cutPvalues[b][i] == -1 )
				{
					fprintf(param->fout,"-");
				}
				else
				{
					fprintf(param->fout,"%f", param->cutPvalues[b][i]);
				}
				if ( b < nlist->size - 1 ) fprintf( param->fout, "," );
			}
			fprintf(param->fout,")\t" );
		}

		if ( param->likelihoodRatios )
		{
			fprintf(param->fout,"(" );
			for ( b = 0 ; b < nlist->size; b++ )
			{
				if ( param->likelihoodRatios[b][i] == -1 )
				{
					fprintf(param->fout,"-");
				}
				else
				{
					fprintf(param->fout,"%f", param->likelihoodRatios[b][i]);
				}
				if ( b < nlist->size - 1 ) fprintf( param->fout, "," );
			}
			fprintf(param->fout,")" );
		}
		fprintf(param->fout, "\n");
	}
}

/**************************************************************************
 * Save.MetaPost
**************************************************************************/


double cafe_report_mp_remark(pString pstr, pTree ptree, pMetapostConfig pmc, va_list ap1)
{
  va_list ap;
  va_copy(ap, ap1);
	double last = cafe_tree_mp_remark(pstr, ptree, pmc, ap);
	//char* title = va_arg(ap,char*);
	va_arg(ap,char*);   // pass first
	pCafeParam param = va_arg(ap, pCafeParam );
	int fid = va_arg(ap,int);
	string_fadd( pstr, "label(btex Max p-value = %4.3f etex, (0.1u, %fu));\n", 
			    param->maximumPvalues[fid],  last - 0.2 );
  va_end(ap);
	return 0;
}

double cafe_report_mp_annotation(pString pstr, pTreeNode pnode, pMetapostConfig pmc, va_list ap1)
{
  va_list ap;
  va_copy(ap, ap1);
	double last = cafe_tree_mp_annotation(pstr, pnode, pmc, ap);
	int idx = pnode->id;

	//char* title = va_arg(ap,char*);
	va_arg(ap,char*); // pass first
	pCafeParam  param = va_arg(ap, pCafeParam );
	int fid = va_arg(ap,int);

	if ( param->cutPvalues && param->cutPvalues[idx][fid] != -1 )
	{
		last -= 0.15;
		string_fadd( pstr, "label.rt( btex bc = %4.3f ", param->cutPvalues[idx][fid] );
		string_fadd( pstr, "etex, mid[%d] + (0,%fu));\n",  pnode->id, last );
	}
	if ( param->likelihoodRatios && param->likelihoodRatios[idx][fid] != -1 )
	{
		last -= 0.15;
		string_fadd( pstr, "label.rt( btex lh = %4.3f ", param->likelihoodRatios[idx][fid]);
		string_fadd( pstr, "etex, mid[%d] + (0,%fu));\n",  pnode->id, last );
	}
  va_end(ap);
	return last;
}


////////////////////////////////////////////////
double cafe_report_mp_summary(pString pstr, pTreeNode pnode, pMetapostConfig pmc, va_list ap1 )
{
  va_list ap;
  va_copy(ap, ap1);
  string_add(pstr, ";\n"); 
  double last = 0;

  //char* title = va_arg(ap,char*);
  va_arg(ap,char*); // pass first
  pCafeParam  param = va_arg(ap, pCafeParam );
  int fid = pnode->id;

  if ( pnode->parent )
  {
    int pid = pnode->parent->id;
    int idx = pid/2 * 2 + (pid - fid > 0 ? 0 : 1); 
    string_fadd( pstr, "xpart mid[%d] = xpart(p[%d]);\n", pnode->id, pnode->id );
    string_fadd( pstr, "ypart mid[%d] = (ypart(p[%d])+ypart(p[%d]))/2;\n", pnode->id, pnode->id, pnode->parent->id );
    last -= 0.15;
    string_fadd( pstr, "label.rt( btex %1.6f ", param->averageExpansion[idx] );
    string_fadd( pstr, "etex, mid[%d] + (0,%fu));\n",  pnode->id, last ); 
    last -= 0.15;
    string_fadd( pstr, "label.rt( btex %4.1f ", ((pCafeNode)pnode)->super.branchlength );    
    string_fadd( pstr, "etex, mid[%d] + (0,%fu));\n",  pnode->id, last );
    last -= 0.15;    
    string_fadd( pstr, "label.rt( btex %f ", ((pCafeNode)pnode)->lambda);
    string_fadd( pstr, "etex, mid[%d] + (0,%fu));\n",  pnode->id, last );
  }
  va_end(ap);
  return last;
}

void cafe_report_summary_metapost(pTree pstree, pCafeParam param,  int id, char* title, double width, double height )
{
  MetapostConfig mc;
  mc.id = id; 
  mc.unit = MP_UNIT_IN;
  mc.dir = MP_DIR_VERTICAL;
  mc.shape = MP_SHAPE_RECT | MP_SHAPE_MOST_CENTER;
  mc.fmod = cafe_report_mp_summary;
  mc.fremark = cafe_report_mp_remark;
  mc.width = width;
  mc.height = height; 
  pString pstr = phylogeny_to_mp( (pTree)pstree, &mc, title, param, param->pfamily->flist->size );
  fprintf( param->fout, "%s\n", pstr->buf );
}
///////////////////////////////////


void cafe_report_pdf(pCafeParam param)
{
	int i;
	pCafeTree pcafe = param->pcafe;
	pArrayList pflist = param->pfamily->flist;
	pString pstr = NULL; 

	pMetapostConfig pmc = cafe_tree_get_default_mpconfig(0, 7, 7);
	pmc->fmod = cafe_report_mp_annotation;
	pmc->fremark = cafe_report_mp_remark;
	for ( i = 0 ; i < pflist->size ; i++ )
	{
		pmc->id = i+1;
		cafe_report_set_viterbi(param,i);
		pCafeFamilyItem pitem = (pCafeFamilyItem)pflist->array[i];
		pstr = phylogeny_to_mp( (pTree)pcafe, pmc, pitem->id, param, i );
		fprintf( param->fout, "%s\n\n", pstr->buf );
		string_free(pstr);
	}
	cafe_report_summary_metapost((pTree)pcafe, param, pflist->size+1, "Summary", 12, 12);
}

void cafe_report(pCafeParam param, int method)
{
	if ( method == CAFE_REPORT_TEXT ) cafe_report_text(param);
	else if ( method == CAFE_REPORT_PDF ) cafe_report_pdf(param);
}


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

int cafe_report_retrieve_data(char* file, pCafeParam param)
{
	int i, j;
	int fs[2] = { 0, 1 };
	int rs[2] = { 1, 2 };

	FILE* fp = fopen(file,"r");		
	if ( fp == NULL )
	{
		fprintf( stderr, "%s: Cannot open file: %s\n", __FUNCTION__, file );
		return -1;
	}
	pString pstr = string_new();
	int nnodes = 0;
	param->expandRemainDecrease = (int**)memory_new(3,sizeof(int*));

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
			param->pcafe = cafe_tree_new(data, fs, rs, 0, 0);
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
			param->averageExpansion = cafe_report_load_data_double_pairs(data, '\t');
		}
		else if ( strncasecmp( pstr->buf, "expansion", 9 ) == 0 )
		{
			param->expandRemainDecrease[0] = cafe_report_load_data_int_pairs(data,'\t');
		}
		else if ( strncasecmp( pstr->buf, "remain", 6 ) == 0 )
		{
			param->expandRemainDecrease[1] = cafe_report_load_data_int_pairs(data,'\t');
		}
		else if ( strncasecmp( pstr->buf, "decrease", 8 ) == 0 )
		{
			param->expandRemainDecrease[2] = cafe_report_load_data_int_pairs(data,'\t');
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

	param->viterbiPvalues = (double**)memory_new_2dim( nnodes-1,num_families, sizeof(double) );
	param->maximumPvalues = (double*)memory_new( num_families, sizeof(double) );
	param->viterbiNodeFamilysizes = (int**)memory_new_2dim( nnodes-1, num_families, sizeof(double) );

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
		param->cutPvalues = (double**)memory_new_2dim(nnodes, num_families, sizeof(double) );		
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
		pitem->pbdc_array = NULL;
		arraylist_add(pcf->flist,pitem);
		sscanf((char*)data->array[2], "%lf", &param->maximumPvalues[i]); 

		pCafeTree ptree = cafe_tree_new( (char*)data->array[1], fs, rs, 0, 0 );
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
			param->viterbiNodeFamilysizes[j/2][i] = pcnode->familysize;
		}

		cafe_report_load_viterbi_pvalue((char*)data->array[3], param->viterbiPvalues, i, nnodes );
		if ( bexist[0] )
		{
			cafe_report_load_bc_or_lhr_list((char*)data->array[4], 
					    param->cutPvalues, i, nnodes );
			param->cutPvalues[param->pcafe->super.root->id][i] = -1.0;
		}
		if ( bexist[1] )
		{
			int idx =  bexist[0] ? 5 : 4;
			cafe_report_load_bc_or_lhr_list((char*)data->array[idx], 
					    param->likelihoodRatios, i, nnodes );
		}
		cafe_tree_free(ptree);
	}
	
	rs[0] = 1; 
	rs[1] = rint(max_size*1.25);
	fs[0] = 0;
    fs[1] =	max_size  + MAX(50,max_size/50);
	cafe_tree_set_parameters(param->pcafe, fs, rs, param->lambda[0] );
	arraylist_free(plines, free);
	return 0;
}
