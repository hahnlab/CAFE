#include "cafe.h"
#include<stdlib.h>
#include<stdio.h>
#include<mathfunc.h>
#include<memalloc.h>
#include<utils.h>

void __cafe_famliy_check_the_pattern(pCafeFamily pcf)
{
	int i, j, k;
	int size = pcf->flist->size;
	for ( i = 0 ; i < size ; i++ )
	{
		pCafeFamilyItem pref = (pCafeFamilyItem)pcf->flist->array[i];
		if ( pref->ref != -1 ) continue;
		pref->ref = i;
		pref->holder = 1;
		for ( j = i + 1 ; j < size ; j++ )
		{
			pCafeFamilyItem pdest = (pCafeFamilyItem)pcf->flist->array[j];
			if ( pdest->ref != -1 ) continue;
			for ( k = 0 ; k < pcf->num_species ; k++ )
			{
				if ( pref->count[k] != pdest->count[k] ) break;
			}
			if ( k == pcf->num_species )
			{
				pdest->ref = i;		
				pdest->holder = 0;
			}
		}
	}
}

int cafe_family_split_cvfiles_byfamily(pCafeParam param, int cv_fold)
{
	int f;
	char* file = param->str_fdata->buf;
	pCafeFamily pcf = param->pfamily;
	int fold = cv_fold;

	for (f=0; f<fold; f++) 
	{
		char trainfile[STRING_BUF_SIZE];
		char validfile[STRING_BUF_SIZE];
		char queryfile[STRING_BUF_SIZE];
		sprintf(trainfile, "%s.%d.train", file, f+1);
		sprintf(validfile, "%s.%d.valid", file, f+1);
		sprintf(queryfile, "%s.%d.query", file, f+1);
		
		FILE* fptrain;
		if ( ( fptrain = fopen( trainfile , "w" ) ) == NULL )
		{
			fprintf(stderr,"ERROR(save): Cannot open %s in write mode.\n", trainfile );
			return -1;
		}
		FILE* fpvalid;
		if ( ( fpvalid = fopen( validfile , "w" ) ) == NULL )
		{
			fprintf(stderr,"ERROR(save): Cannot open %s in write mode.\n", trainfile );
			return -1;
		}
		FILE* fpquery;
		if ( ( fpquery = fopen( queryfile , "w" ) ) == NULL )
		{
			fprintf(stderr,"ERROR(save): Cannot open %s in write mode.\n", queryfile );
			return -1;
		}
		
		int i, n;
		char buf[STRING_STEP_SIZE];
		string_pchar_join(buf,"\t", pcf->num_species, pcf->species );
		fprintf(fptrain,"Desc\tFamily ID\t%s\n", buf );	
		fprintf(fpvalid,"Desc\tFamily ID\t%s\n", buf );	
		fprintf(fpquery,"Desc\tFamily ID\tspecies:count\n");	
		
		for ( i = 0 ; i < pcf->flist->size ; i++ )
		{
			FILE* fp;
			pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[i];
			if ((f+i+1)%fold == 0) {
				// randomly pick one leaf and print for query
				int randomidx = (int)(unifrnd()*pcf->num_species);
				fprintf(fpquery,"%s\t%s\t%s\t%d\n", pitem->desc, pitem->id, pcf->species[randomidx], pitem->count[randomidx]);	
				//
				fp = fpvalid;
			}
			else {
				fp = fptrain;
			}
			fprintf(fp,"%s\t%s\t%d", pitem->desc,pitem->id, pitem->count[0]);	
			for ( n =  1 ; n < pcf->num_species ; n++ )
			{
				fprintf(fp,"\t%d", pitem->count[n]);	
			}
			fprintf(fp,"\n");
			
		}
		fclose(fptrain);
		fclose(fpvalid);
		fclose(fpquery);
	}
	return 0;
}

void cafe_family_clean_cvfiles_byfamily(pCafeParam param, int cv_fold) 
{
	int i;
	char buf[STRING_BUF_SIZE];
	char* file = param->str_fdata->buf;
	FILE* fp = fopen(param->str_fdata->buf,"r");
	if ( fp == NULL )
	{
		fprintf( stderr, "Cannot open family file: %s\n", file );
		//		print_error(__FILE__,(char*)__FUNCTION__,__LINE__, "Cannot open file: %s ", file );
		return;
	}
	if ( fgets(buf,STRING_BUF_SIZE,fp) == NULL )
	{
		fclose(fp);
		fprintf( stderr, "Empty family file: %s\n", file );
		//		print_error(__FILE__,(char*)__FUNCTION__,__LINE__, "Input file is empty" );
		return;
	}
	string_pchar_chomp(buf);
	pArrayList data = string_pchar_split( buf, '\t');
	fclose(fp);
	
	for (i=0; i<cv_fold; i++) {
		char tmp[STRING_BUF_SIZE];
		sprintf(tmp, "%s.%d.train", file, i+1);
		if (remove(tmp) == -1) { 
			fprintf(stderr, "Error in deleting file %s", tmp);
			return;
		}
		sprintf(tmp, "%s.%d.valid", file, i+1);
		if (remove(tmp) == -1) { 
			fprintf(stderr, "Error in deleting file %s", tmp);
			return;
		}		
		sprintf(tmp, "%s.%d.query", file, i+1);
		if (remove(tmp) == -1) { 
			fprintf(stderr, "Error in deleting file %s", tmp);
			return;
		}		
	}
	arraylist_free(data, free);	
	// free cv_test_count_list and cv_species_count_list
	arraylist_free(param->cv_test_count_list, free);
	arraylist_free(param->cv_test_species_list, free);
	
}


void cafe_family_split_cvfiles_byspecies(pCafeParam param)
{
	int i, j;
	char buf[STRING_BUF_SIZE];
	char* file = param->str_fdata->buf;
	FILE* fp = fopen(param->str_fdata->buf,"r");
	if ( fp == NULL )
	{
		fprintf( stderr, "Cannot open family file: %s\n", file );
		//		print_error(__FILE__,(char*)__FUNCTION__,__LINE__, "Cannot open file: %s ", file );
		return;
	}
	if ( fgets(buf,STRING_BUF_SIZE,fp) == NULL )
	{
		fclose(fp);
		fprintf( stderr, "Empty family file: %s\n", file );
		//		print_error(__FILE__,(char*)__FUNCTION__,__LINE__, "Input file is empty" );
		return;
	}
	string_pchar_chomp(buf);
	pArrayList data = string_pchar_split( buf, '\t');
	int species_num = data->size-2;
	fclose(fp);
	
	for (j=species_num+2; j>2; j--) {
		pString command = string_new_with_string( "awk '{print $1\"\\t\"$2" );
		for (i=3; i<3+species_num; i++) {
			if (i!=j) {
				char tmp[STRING_BUF_SIZE];
				sprintf(tmp,"\"\\t\"$%d", i);
				string_add(command, tmp);
			}
		}
		char tmp[STRING_BUF_SIZE];
		sprintf(tmp, "}' %s > %s.%s.train", file, file, (char*)data->array[j-1]);
		string_add(command, tmp);
		//printf("%s\n", command->buf);
		system(command->buf);
		string_free(command);
		pString command2 = string_new_with_string( "awk '{print $1\"\\t\"$2" );
		sprintf(tmp, "\"\\t\"$%d", j);
		string_add(command2, tmp);
		sprintf(tmp, "}' %s > %s.%s.valid", file, file, (char*)data->array[j-1]);
		string_add(command2, tmp);
		//printf("%s\n", command2->buf);
		system(command2->buf);
		string_free(command2);
		
	}
	arraylist_free(data, free);
}


void cafe_family_clean_cvfiles_byspecies(pCafeParam param) 
{
	int j;
	char buf[STRING_BUF_SIZE];
	char* file = param->str_fdata->buf;
	FILE* fp = fopen(param->str_fdata->buf,"r");
	if ( fp == NULL )
	{
		fprintf( stderr, "Cannot open family file: %s\n", file );
		//		print_error(__FILE__,(char*)__FUNCTION__,__LINE__, "Cannot open file: %s ", file );
		return;
	}
	if ( fgets(buf,STRING_BUF_SIZE,fp) == NULL )
	{
		fclose(fp);
		fprintf( stderr, "Empty family file: %s\n", file );
		//		print_error(__FILE__,(char*)__FUNCTION__,__LINE__, "Input file is empty" );
		return;
	}
	string_pchar_chomp(buf);
	pArrayList data = string_pchar_split( buf, '\t');
	int species_num = data->size-2;
	fclose(fp);
	
	for (j=species_num+2; j>2; j--) {
		char tmp[STRING_BUF_SIZE];
		sprintf(tmp, "%s.%s.train", file, (char*)data->array[j-1]);
		if (remove(tmp) == -1) { 
			fprintf(stderr, "Error in deleting file %s", tmp);
			return;
		}
		sprintf(tmp, "%s.%s.valid", file, (char*)data->array[j-1]);
		if (remove(tmp) == -1) { 
			fprintf(stderr, "Error in deleting file %s", tmp);
			return;
		}		
	}
	arraylist_free(data, free);	
	// free cv_test_count_list
	arraylist_free(param->cv_test_count_list, free);
	
}

void cafe_family_read_query_family(pCafeParam param, const char* file)
{
	FILE* fp = fopen(file,"r");
	char buf[STRING_BUF_SIZE];
	if ( fp == NULL )
	{
		fprintf( stderr, "Cannot open file: %s\n", file );
		//		print_error(__FILE__,(char*)__FUNCTION__,__LINE__, "Cannot open file: %s ", file );
		return;
	}
	int i;
	if ( fgets(buf,STRING_BUF_SIZE,fp) == NULL )
	{
		fclose(fp);
		fprintf( stderr, "Empty file: %s\n", file );
		//		print_error(__FILE__,(char*)__FUNCTION__,__LINE__, "Input file is empty" );
		return;
	}
	string_pchar_chomp(buf);
	pArrayList data = string_pchar_split( buf, '\t');
	param->cv_test_species_list = arraylist_new(11000);
	param->cv_test_count_list = arraylist_new(11000);
	
	int maxsize = 0;
	for(  i = 0 ; fgets(buf,STRING_BUF_SIZE,fp) ; i++ )	
	{
		char* species = NULL;
		int* count = memory_new(1, sizeof(int));
		//		string_pchar_chomp(buf);
		pArrayList data = string_pchar_split( buf, '\t');
		if ( data->size != 4 )
		{
			fprintf(stderr,"format needs to be 3 columns ( desc id familysize ): %s\n",(char*)data->array[i] );
		}
		//// this part is reading the counts
		if ( data->array[2] && data->array[3] )
		{
			species = strdup(data->array[2]);
			*count = atoi((char*)data->array[3]);
			if ( maxsize < *count )
			{
				maxsize = *count;
			}
		}
		else
		{
			*count = 0;
		}
		/// end read count
		arraylist_free(data,free);
		arraylist_add(param->cv_test_species_list,(void*)species);
		arraylist_add(param->cv_test_count_list,(void*)count);
	}
	arraylist_free(data,free);
	fclose(fp);
	return;
}



void cafe_family_read_validate_species(pCafeParam param, const char* file)
{
	FILE* fp = fopen(file,"r");
	char buf[STRING_BUF_SIZE];
	if ( fp == NULL )
	{
		fprintf( stderr, "Cannot open file: %s\n", file );
		//		print_error(__FILE__,(char*)__FUNCTION__,__LINE__, "Cannot open file: %s ", file );
		return;
	}
	int i;
	if ( fgets(buf,STRING_BUF_SIZE,fp) == NULL )
	{
		fclose(fp);
		fprintf( stderr, "Empty file: %s\n", file );
		//		print_error(__FILE__,(char*)__FUNCTION__,__LINE__, "Input file is empty" );
		return;
	}
	string_pchar_chomp(buf);
	pArrayList data = string_pchar_split( buf, '\t');
	param->cv_species_name = strdup(data->array[2]);
	param->cv_test_count_list = arraylist_new(11000);
	
	int maxsize = 0;
	for(  i = 0 ; fgets(buf,STRING_BUF_SIZE,fp) ; i++ )	
	{
		int* count = memory_new(1, sizeof(int));
		//		string_pchar_chomp(buf);
		pArrayList data = string_pchar_split( buf, '\t');
		if ( data->size != 3 )
		{
			fprintf(stderr,"format needs to be 3 columns ( desc id familysize ): %s\n",(char*)data->array[i] );
		}
		//// this part is reading the counts
		if ( data->array[2] )
		{
			*count = atoi((char*)data->array[2]);
			if ( maxsize < *count )
			{
				maxsize = *count;
			}
		}
		else
		{
			*count = 0;
		}
		/// end read count
		arraylist_free(data,free);
		arraylist_add(param->cv_test_count_list,(void*)count);
	}
	arraylist_free(data,free);
	fclose(fp);
	return;
}

/*! \brief Synchronize a family and tree together
*
*  Sets the value of index in pcf to the node index in the tree of the 
*  leaf with the matching species ID. 
*  \returns 0 on success, or -1 if there is a species in the tree with
*  no matching species in the family
*/int cafe_family_set_species_index(pCafeFamily pcf, pCafeTree pcafe )
{
	int i,j;
	pTree ptree = (pTree)pcafe;
	for ( i = 0 ; i < pcf->num_species ;  i++ )
	{
		if (pcf->species[i][0] == '-') {
			//pPhylogenyNode pnode = (pPhylogenyNode)ptree->nlist->array[j];
			int nodeidx = atoi(&pcf->species[i][1]);
			pcf->index[i] = nodeidx;
		}
		else {
			for( j = 0 ; j < ptree->nlist->size; j+=2 )		// find node with same name among leaves
			{
				pPhylogenyNode pnode = (pPhylogenyNode)ptree->nlist->array[j];
				if ( pnode->name[0] & 0x80 ) {
					continue;	// meta-bit set: skip already processed 
				}
				if ( string_pchar_cmp_ignore_case(pnode->name, pcf->species[i] ) )
				{
					pnode->name[0] |= 0x80;	// set meta-bit.
					pcf->index[i] = j;
					break;
				}
			}
		}
		//printf("%d: %s => %d\n", i, pcf->species[i], pcf->index[i] );
	}

	int err = 0;
	for( j = 0 ; j < ptree->nlist->size; j+=2 )
	{
		pPhylogenyNode pnode = (pPhylogenyNode)ptree->nlist->array[j];
		if ( !(pnode->name[0] & 0x80 ) )
		{
			fprintf(stderr,"Family has no %s species\n", pnode->name );
			err = -1;
		}
	}

	for( j = 0 ; j < ptree->nlist->size; j+=2 )
	{
		pPhylogenyNode pnode = (pPhylogenyNode)ptree->nlist->array[j];
		pnode->name[0] &= 0x7F;
	}
	return err;
}


int cafe_family_get_species_index(pCafeFamily pcf, char* speciesname) 
{
    int i;
    int returnidx = -1;
	for ( i = 0 ; i < pcf->num_species ;  i++ )
	{
        if ( string_pchar_cmp_ignore_case(speciesname, pcf->species[i] ) )
        {
            returnidx = pcf->index[i];
            break;
        }
    }
    return returnidx;
}

/*! \brief Copy sizes of an individual gene family into the tree
*
* Takes the given gene family and copies the size of the family in the species
* into the tree node that corresponds to that species
*
* Side effect is to set the familysize member for each node in the tree, but
* since that is a convenience variable it shouldn't make any difference
*/
void cafe_family_set_size(pCafeFamily pcf, int idx, pCafeTree pcafe)
{
	assert(idx < pcf->flist->size);

	pTree ptree = (pTree)pcafe;
	for (int i = 0; i< ptree->nlist->size; i++) {
		((pCafeNode)ptree->nlist->array[i])->familysize = -1;
	}
	pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[idx];
	for (int i = 0 ; i < pcf->num_species ; i++ )
	{
		if (pcf->index[i] < 0)
		{
			fprintf(stderr, "Warning: Tree and family indices not synchronized\n");
			continue;
		}
		if (pcf->index[i] >= ptree->size)
		{
			fprintf(stderr, "Inconsistency in tree size");
			exit(-1);
		}
		pCafeNode pcnode = (pCafeNode)ptree->nlist->array[pcf->index[i]];
		pcnode->familysize = pitem->count[i];
	}
}

void cafe_family_set_size_with_family_forced(pCafeFamily pcf, int idx, pCafeTree pcafe)
{
	cafe_family_set_size(pcf, idx, pcafe);
	pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[idx];
	int i;
	int max = 0;
	for ( i = 0 ; i < pcf->num_species ; i++ )
	{
		if( pcf->index[i] < 0 ) continue;
		if ( max < pitem->count[i] )
		{
			max = pitem->count[i];
		}
	}
	pcafe->rootfamilysizes[0] = 1;
	pcafe->rootfamilysizes[1] = rint(max * 1.25);
	pcafe->familysizes[1] = max + MAX(50,max/5);
	pcafe->rfsize = pcafe->rootfamilysizes[1] - pcafe->rootfamilysizes[0] + 1;
}

void cafe_family_set_size_with_family(pCafeFamily pcf, int idx, pCafeTree pcafe )
{
	cafe_family_set_size(pcf, idx, pcafe);			// set the size of leaves according to data
	pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[idx];
	int i;
	int max = 0;
	for ( i = 0 ; i < pcf->num_species ; i++ )
	{
		if( pcf->index[i] < 0 ) continue;
		if ( max < pitem->count[i] )
		{
			max = pitem->count[i];
		}
	}
	
	if ( pitem->maxlh >= 0 )			// adjusting the family size based on the maxlh index
	{
		//pcafe->rootfamilysizes[0] = MAX(1,pitem->maxlh - 10);
		pcafe->rootfamilysizes[0] = 1;
		pcafe->rootfamilysizes[1] = pitem->maxlh + 20;
		pcafe->familysizes[1] = MAX( MAX(pitem->maxlh + 50, rint(pitem->maxlh * 1.2)), 
					                 MAX(max + 50, rint(max*1.2) ) );
		pcafe->rfsize = pcafe->rootfamilysizes[1] - pcafe->rootfamilysizes[0] + 1;
//		max = MAX( pcafe->rootfamilysizes[1], pcafe->familysizes[1] );
	}
	
}

void cafe_family_set_truesize_with_family(pCafeFamily pcf, int idx, pCafeTree pcafe )
{
	int i;
	int sumoffamilysize = 0;
	// set the size of leaves according to data
	pTree ptree = (pTree)pcafe;
	for ( i = 0; i< ptree->nlist->size; i++) {
		((pCafeNode)ptree->nlist->array[i])->familysize = -1;
	}
	pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[idx];
	for ( i = 0 ; i < pcf->num_species ; i++ )
	{
		if( pcf->index[i] < 0 ) continue;
		pCafeNode pcnode = (pCafeNode)ptree->nlist->array[pcf->index[i]];
		pcnode->familysize = pitem->count[i];
		sumoffamilysize += pitem->count[i];
	}
	if (sumoffamilysize == 0) 
	{
		//fprintf(stderr, "empty family\n");
	}
	// adjust family size range
	int max = 0;
	for ( i = 0 ; i < pcf->num_species ; i++ )
	{
		if( pcf->index[i] < 0 ) continue;
		if ( max < pitem->count[i] )
		{
			max = pitem->count[i];
		}
	}	
	if ( pitem->maxlh >= 0 )			// adjusting the family size based on the maxlh index
	{
		//pcafe->rootfamilysizes[0] = MAX(1,pitem->maxlh - 10);
		pcafe->rootfamilysizes[0] = 1;
		pcafe->rootfamilysizes[1] = pitem->maxlh + 20;
		pcafe->familysizes[1] = MAX( MAX(pitem->maxlh + 50, rint(pitem->maxlh * 1.2)), 
									MAX(max + 50, rint(max*1.2) ) );
		pcafe->rfsize = pcafe->rootfamilysizes[1] - pcafe->rootfamilysizes[0] + 1;
		//		max = MAX( pcafe->rootfamilysizes[1], pcafe->familysizes[1] );
	}
	
}



void cafe_family_set_size_by_species(char* speciesname, int size, pCafeTree pcafe)
{
	pTree ptree = (pTree)pcafe;
	int i;
	// only look for leaf nodes
	for ( i = 0; i< ptree->nlist->size; i++) {
		((pCafeNode)ptree->nlist->array[i])->familysize = -1;
	}
	for ( i = 0; i< ptree->nlist->size; i=i+2) {
		if (strcmp(((pPhylogenyNode)ptree->nlist->array[i])->name, speciesname) == 0) {
			((pCafeNode)ptree->nlist->array[i])->familysize = size;
		}
	}
}



void cafe_family_reset_maxlh(pCafeFamily pcf)
{
	int i;
	pArrayList flist = pcf->flist;
	for ( i = 0 ; i < flist->size ; i++ )
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)flist->array[i];			
		pitem->maxlh = -1;
	}
}

int cafe_family_get_index(pCafeFamily pcf, const char* szid)
{
	int i;
	pArrayList flist = pcf->flist;
	for ( i = 0 ; i < flist->size ; i++ )
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)flist->array[i];			
		if ( pitem->id && strcmp( pitem->id, szid ) == 0 ) break;
	}
	return flist->size == i ? -1 : i;
}

pCafeFamilyItem cafe_family_get_family_item(pCafeFamily pcf, const char* szid )
{
	int i = cafe_family_get_index( pcf, szid );	
	if ( i == -1 ) return NULL;
	return (pCafeFamilyItem)pcf->flist->array[i];
}

/*
                 p      r
      0   0   => 0      0 
	  0   1   => 1/0    0 ( 0 side has child with 1, then 1 )
      0   1/0 => 0      0
      1   1/0 => 1      1
      1   1   => 1      1
      1/0 1/0 => 1/0    1
 */

void init_family_size(family_size_range* fs, int max)
{
	fs->root_min = 1;	// Must be 1, not 0
	fs->root_max = MAX(30, rint(max * 1.25));

	fs->max = max + MAX(50, max / 5);
	fs->min = 0;
}






