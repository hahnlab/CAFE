#include <string>
#include <stdexcept>
#include <fstream>

#include "cross_validator.h"
#include "gene_family.h"

void cross_validator::clean_by_family(std::string file, int cv_fold)
{
  int i;
//  char* file = param->str_fdata->buf;
  for (i = 0; i<cv_fold; i++) {
    char tmp[STRING_BUF_SIZE];
    sprintf(tmp, "%s.%d.train", file.c_str(), i + 1);
    if (remove(tmp) == -1) {
      fprintf(stderr, "Error in deleting file %s", tmp);
      return;
    }
    sprintf(tmp, "%s.%d.valid", file.c_str(), i + 1);
    if (remove(tmp) == -1) {
      fprintf(stderr, "Error in deleting file %s", tmp);
      return;
    }
    sprintf(tmp, "%s.%d.query", file.c_str(), i + 1);
    if (remove(tmp) == -1) {
      fprintf(stderr, "Error in deleting file %s", tmp);
      return;
    }
  }
  // free cv_test_count_list and cv_species_count_list
  arraylist_free(cv_test_count_list, free);
  arraylist_free(cv_test_species_list, free);

}


void cross_validator::clean_by_species(std::string file)
{
  int j;
  char buf[STRING_BUF_SIZE];
  FILE* fp = fopen(file.c_str(), "r");
  if (fp == NULL)
  {
    fprintf(stderr, "Cannot open family file: %s\n", file.c_str());
    //		print_error(__FILE__,(char*)__FUNCTION__,__LINE__, "Cannot open file: %s ", file );
    return;
  }
  if (fgets(buf, STRING_BUF_SIZE, fp) == NULL)
  {
    fclose(fp);
    fprintf(stderr, "Empty family file: %s\n", file.c_str());
    //		print_error(__FILE__,(char*)__FUNCTION__,__LINE__, "Input file is empty" );
    return;
  }
  string_pchar_chomp(buf);
  pArrayList data = string_pchar_split(buf, '\t');
  int species_num = data->size - 2;
  fclose(fp);

  for (j = species_num + 2; j>2; j--) {
    char tmp[STRING_BUF_SIZE];
    sprintf(tmp, "%s.%s.train", file.c_str(), (char*)data->array[j - 1]);
    if (remove(tmp) == -1) {
      fprintf(stderr, "Error in deleting file %s", tmp);
      return;
    }
    sprintf(tmp, "%s.%s.valid", file.c_str(), (char*)data->array[j - 1]);
    if (remove(tmp) == -1) {
      fprintf(stderr, "Error in deleting file %s", tmp);
      return;
    }
  }
  arraylist_free(data, free);
  // free cv_test_count_list
  arraylist_free(cv_test_count_list, free);

}

void cross_validator::read_validate_species(const char* file)
{
  FILE* fp = fopen(file, "r");
  char buf[STRING_BUF_SIZE];
  if (fp == NULL)
  {
    fprintf(stderr, "Cannot open file: %s\n", file);
    //		print_error(__FILE__,(char*)__FUNCTION__,__LINE__, "Cannot open file: %s ", file );
    return;
  }
  int i;
  if (fgets(buf, STRING_BUF_SIZE, fp) == NULL)
  {
    fclose(fp);
    fprintf(stderr, "Empty file: %s\n", file);
    //		print_error(__FILE__,(char*)__FUNCTION__,__LINE__, "Input file is empty" );
    return;
  }
  string_pchar_chomp(buf);
  pArrayList data = string_pchar_split(buf, '\t');
  cv_species_name = strdup((char *)data->array[2]);
  cv_test_count_list = arraylist_new(11000);

  int maxsize = 0;
  for (i = 0; fgets(buf, STRING_BUF_SIZE, fp); i++)
  {
    int* count = (int *)memory_new(1, sizeof(int));
    //		string_pchar_chomp(buf);
    pArrayList data = string_pchar_split(buf, '\t');
    if (data->size != 3)
    {
      fprintf(stderr, "format needs to be 3 columns ( desc id familysize ): %s\n", (char*)data->array[i]);
    }
    //// this part is reading the counts
    if (data->array[2])
    {
      *count = atoi((char*)data->array[2]);
      if (maxsize < *count)
      {
        maxsize = *count;
      }
    }
    else
    {
      *count = 0;
    }
    /// end read count
    arraylist_free(data, free);
    arraylist_add(cv_test_count_list, (void*)count);
  }
  arraylist_free(data, free);
  fclose(fp);
  return;
}

void cross_validator::read_query_family(const char* file)
{
  FILE* fp = fopen(file, "r");
  char buf[STRING_BUF_SIZE];
  if (fp == NULL)
  {
    fprintf(stderr, "Cannot open file: %s\n", file);
    //		print_error(__FILE__,(char*)__FUNCTION__,__LINE__, "Cannot open file: %s ", file );
    return;
  }
  int i;
  if (fgets(buf, STRING_BUF_SIZE, fp) == NULL)
  {
    fclose(fp);
    fprintf(stderr, "Empty file: %s\n", file);
    //		print_error(__FILE__,(char*)__FUNCTION__,__LINE__, "Input file is empty" );
    return;
  }
  string_pchar_chomp(buf);
  pArrayList data = string_pchar_split(buf, '\t');
  cv_test_species_list = arraylist_new(11000);
  cv_test_count_list = arraylist_new(11000);

  int maxsize = 0;
  for (i = 0; fgets(buf, STRING_BUF_SIZE, fp); i++)
  {
    char* species = NULL;
    int* count = (int *)memory_new(1, sizeof(int));
    //		string_pchar_chomp(buf);
    pArrayList data = string_pchar_split(buf, '\t');
    if (data->size != 4)
    {
      fprintf(stderr, "format needs to be 3 columns ( desc id familysize ): %s\n", (char*)data->array[i]);
    }
    //// this part is reading the counts
    if (data->array[2] && data->array[3])
    {
      species = strdup((char *)data->array[2]);
      *count = atoi((char*)data->array[3]);
      if (maxsize < *count)
      {
        maxsize = *count;
      }
    }
    else
    {
      *count = 0;
    }
    /// end read count
    arraylist_free(data, free);
    arraylist_add(cv_test_species_list, (void*)species);
    arraylist_add(cv_test_count_list, (void*)count);
  }
  arraylist_free(data, free);
  fclose(fp);
  return;
}

double cross_validator::validate_by_species(pCafeParam param, const char* validatefile, const char* errortype)
{
  int i, j;
  read_validate_species(validatefile);
  if (cv_test_count_list == NULL) return -1;
  // now compare reconstructed count to true count	
  pCafeTree pcafe = param->pcafe;

  reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, &param->family_size);
  pArrayList estimate_size = arraylist_new(cv_test_count_list->size);
  for (i = 0; i< param->pfamily->flist->size; i++)
  {
    pCafeFamilyItem pitem = (pCafeFamilyItem)param->pfamily->flist->array[i];
    cafe_family_set_size(param->pfamily, pitem, pcafe);
    if (param->posterior) {
      cafe_tree_viterbi_posterior(pcafe, param);
    }
    else {
      cafe_tree_viterbi(pcafe);
    }
    for (j = 0; j<pcafe->super.nlist->size; j++) {
      char* nodename = ((pPhylogenyNode)pcafe->super.nlist->array[j])->name;
      if (nodename && (strcmp(nodename, cv_species_name) == 0)) {
        int* pFamilysize = (int *)memory_new(1, sizeof(int));
        *pFamilysize = ((pCafeNode)pcafe->super.nlist->array[j])->familysize;
        arraylist_add(estimate_size, (void*)pFamilysize);
      }

    }
  }
  cafe_free_birthdeath_cache(pcafe);
  double MSE = 0;
  double MAE = 0;
  if (cv_test_count_list->size != param->pfamily->flist->size)
    throw std::runtime_error("list sizes don't match\n");

  for (i = 0; i<cv_test_count_list->size; i++) {
    int error = (*((int*)cv_test_count_list->array[i]) - *((int*)estimate_size->array[i]));
    MSE += pow(error, 2);
    MAE += abs(error);
  }
  MSE = MSE / (cv_test_count_list->size);
  MAE = MAE / (cv_test_count_list->size);
  cafe_log(param, "MSE %f\n", MSE);
  cafe_log(param, "MAE %f\n", MAE);

  arraylist_free(estimate_size, free);

  double returnerror = -1;
  if (strncmp(errortype, "MSE", 3) == 0) {
    returnerror = MSE;
  }
  else if (strncmp(errortype, "MAE", 3) == 0) {
    returnerror = MAE;
  }
  return returnerror;
}

double cross_validator::validate_by_family(pCafeParam param, const char* queryfile, const char* truthfile, const char* errortype)
{
  int i, j;
  double MSE = 0;
  double MAE = 0;
  double SSE = 0;
  double SAE = 0;
  read_query_family(queryfile);
  if (cv_test_count_list == NULL) return -1;

  // read in validation data
  std::ifstream ifst(truthfile);
  pCafeFamily truthfamily = load_gene_families(ifst, 1, '\t');
  if (truthfamily == NULL) {
    fprintf(stderr, "failed to read in true values %s\n", truthfile);
    return -1;
  }

  // now compare reconstructed count to true count	
  pCafeTree pcafe = param->pcafe;
  pCafeTree truthtree = cafe_tree_copy(pcafe);
  // set parameters
  if (truthtree)
  {
    cafe_family_set_species_index(truthfamily, truthtree);
  }

  reset_birthdeath_cache(param->pcafe, param->parameterized_k_value, &param->family_size);

  for (i = 0; i< cv_test_count_list->size; i++)
  {
    int* testcnt = (int*)cv_test_count_list->array[i];
    pCafeFamilyItem pitem = (pCafeFamilyItem)truthfamily->flist->array[i];
    cafe_family_set_size(truthfamily, pitem, truthtree);
    cafe_family_set_size_by_species((char *)cv_test_species_list->array[i], *testcnt, pcafe);
    if (param->posterior) {
      cafe_tree_viterbi_posterior(pcafe, param);
    }
    else {
      cafe_tree_viterbi(pcafe);
    }
    // leaf nodes SSE
    SSE = 0;
    SAE = 0;
    int nodecnt = 0;
    for (j = 0; j<pcafe->super.nlist->size; j = j + 2) {
      int error = ((pCafeNode)truthtree->super.nlist->array[j])->familysize - ((pCafeNode)pcafe->super.nlist->array[j])->familysize;
      SSE += pow(error, 2);
      SAE += abs(error);
      nodecnt++;
    }
    MSE += SSE / nodecnt;
    MSE += SAE / nodecnt;
  }
  cafe_free_birthdeath_cache(pcafe);

  MSE = MSE / cv_test_count_list->size;
  MAE = MAE / cv_test_count_list->size;
  cafe_log(param, "MSE %f\n", MSE);
  cafe_log(param, "MAE %f\n", MSE);

  double returnerror = -1;
  if (strncmp(errortype, "MSE", 3) == 0) {
    returnerror = MSE;
  }
  else if (strncmp(errortype, "MAE", 3) == 0) {
    returnerror = MAE;
  }
  return returnerror;
}


