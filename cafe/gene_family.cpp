#include "gene_family.h"
#include <sstream>

extern "C" {
#include "cafe.h"
extern pCafeParam cafe_param;
void __cafe_famliy_check_the_pattern(pCafeFamily pcf);
}

using namespace std;

vector<string> tokenize(string s)
{
  vector<string> result;
  istringstream iss(s);

  while (iss.good()) {
    string tmp;
    iss >> tmp;
    if (tmp.size() > 0)
      result.push_back(tmp);
  }

  return result;
}


pCafeFamily cafe_family_init(const std::vector<std::string>& species_list)
{
  pCafeFamily pcf = (pCafeFamily)memory_new(1, sizeof(CafeFamily));
  pcf->num_species = species_list.size();
  pcf->max_size = 0;
  pcf->flist = arraylist_new(11000);

  pcf->species = (char**)memory_new(species_list.size(), sizeof(char*));
  for (size_t i = 0; i < species_list.size(); i++)
  {
    pcf->species[i] = new char[species_list[i].size()+1];
    strcpy(pcf->species[i], species_list[i].c_str());
  }

  pcf->index = (int*)memory_new(pcf->num_species, sizeof(int));
  for (int i = 0; i < pcf->num_species; ++i)
  {
    pcf->index[i] = -1;
  }

  return pcf;
}

void cafe_family_free(pCafeFamily pcf)
{
  for (int i = 0; i < pcf->num_species; i++) 
  { 
    delete [] pcf->species[i];  
    pcf->species[i] = NULL;
  }
  memory_free(pcf->species);
  pcf->species = NULL;
  memory_free(pcf->index);
  pcf->index = NULL;
  arraylist_free(pcf->flist, (freefunc)cafe_family_item_free);
  memory_free(pcf);
  pcf = NULL;
}

pCafeFamily cafe_family_new(const char* file, int bpatcheck)
{
  FILE* fp = fopen(file, "r");
  char buf[STRING_BUF_SIZE];
  if (fp == NULL)
  {
    fprintf(stderr, "Cannot open file: %s\n", file);
    return NULL;
  }
  if (fgets(buf, STRING_BUF_SIZE, fp) == NULL)
  {
    fclose(fp);
    fprintf(stderr, "Empty file: %s\n", file);
    return NULL;
  }
  string_pchar_chomp(buf);

  vector<string> species_list = tokenize(buf);
  species_list.erase(species_list.begin(), species_list.begin()+2); // first two items are description and ID - delete them
  pCafeFamily pcf = cafe_family_init(species_list);

  for (int i = 0; fgets(buf, STRING_BUF_SIZE, fp); i++)
  {
    pArrayList data = string_pchar_split(buf, '\t');
    cafe_family_add_item(pcf, data);
    arraylist_free(data, NULL);
  }
  if (bpatcheck)
  {
    __cafe_famliy_check_the_pattern(pcf);
  }
  fclose(fp);
  return pcf;
}


double cross_validate_by_family(const char* queryfile, const char* truthfile, const char* errortype)
{
  int i, j;
  double MSE = 0;
  double MAE = 0;
  double SSE = 0;
  double SAE = 0;
  cafe_family_read_query_family(cafe_param, queryfile);
  if (cafe_param->cv_test_count_list == NULL) return -1;

  // read in validation data
  pCafeFamily truthfamily = cafe_family_new(truthfile, 1);
  if (truthfamily == NULL) {
    fprintf(stderr, "failed to read in true values %s\n", truthfile);
    return -1;
  }

  // now compare reconstructed count to true count	
  pCafeTree pcafe = cafe_param->pcafe;
  pCafeTree truthtree = cafe_tree_copy(pcafe);
  // set parameters
  if (truthtree)
  {
    cafe_family_set_species_index(truthfamily, truthtree);
  }

  reset_birthdeath_cache(cafe_param->pcafe, cafe_param->parameterized_k_value, &cafe_param->family_size);

  for (i = 0; i< cafe_param->cv_test_count_list->size; i++)
  {
    int* testcnt = (int*)cafe_param->cv_test_count_list->array[i];
    cafe_family_set_size(truthfamily, i, truthtree);
    cafe_family_set_size_by_species((char *)cafe_param->cv_test_species_list->array[i], *testcnt, pcafe);
    if (cafe_param->posterior) {
      cafe_tree_viterbi_posterior(pcafe, cafe_param);
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

  MSE = MSE / cafe_param->cv_test_count_list->size;
  MAE = MAE / cafe_param->cv_test_count_list->size;
  cafe_log(cafe_param, "MSE %f\n", MSE);
  cafe_log(cafe_param, "MAE %f\n", MSE);

  double returnerror = -1;
  if (strncmp(errortype, "MSE", 3) == 0) {
    returnerror = MSE;
  }
  else if (strncmp(errortype, "MAE", 3) == 0) {
    returnerror = MAE;
  }
  return returnerror;
}
