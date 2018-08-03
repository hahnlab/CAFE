/*! \page Load Load
* \code{.sh}
# load -i filename [-t # of CPU threads] [-l filename] [-p &alpha;] [-r # of random samples
* \endcode

-i filename: Enter the path to the file containing gene family data. The data file format must be tab-delimited with 
UNIX line endings. Family description may contain spaces (but not tabs). The first line must contain labels in the order: 
* Description
* ID, 
* A list of tab-delimited taxon names

\note If you do not have a Description or ID, CAFE still requires two tabs at the beginning of each line). 

The taxon names must be spelled exactly as they are in the provided phylogenetic tree. Each subsequent line then 
corresponds to a single gene family. If the data file contains taxa that do not appear in the tree structure, 
they are not considered in the analysis.

Here is an example of an input data file:

\code
Description ID Chimp Human Mouse Rat Dog 
EF 1 ALPHA ENSF00000000004 5 8 6 12 40 
HLA CLASS II ENSF00000000007 4 4 3 3 3 
HLA CLASS I ENSF00000000014 5 3 5 6 3 
RAG 1 ENSF00000000015 1 1 1 1 1 
IG HEAVY CHAIN ENSF00000000020 32 42 51 60 18 
ACTIN ENSF00000000027 27 30 22 28 25 
OPSIN ENSF00000000029 2 2 2 2 2 
HEAVY CHAIN ENSF00000000030 25 25 23 24 18
\endcode

If the file is loaded correctly, CAFE will output summary information about the current data file to the log file.

\b -t <em># of CPU threads</em>: The maximum number (integer) of CPU threads to be used. Default: 8.

\b -l <em>filename</em>: Enter the path to the file where CAFE will write the main output. This file will contain a 
summary of input parameters as well as details of &lambda; searches, including likelihood scores and maximum likelihood 
values of &lambda;. If the file does not exist, CAFE will create it for you; if the file already exists, CAFE will append 
the results to the previous file. Default: output to screen (no log file created).

\b -p &alpha;: For each family in the data file, CAFE computes a probability (p-value) of observing the data given the 
average rate of gain and loss of genes. All else being equal, families with more variance in size are expected to have 
lower p-values. The significance level (&alpha;, a float) allows the user to specify the cutoff for subsequent analyses. 
Families with p-values larger than the designated significance level will not be included in the identification of the 
most unlikely branch. Default: 0.01.

\b -r <em># of random samples</em>: To determine the probability of a gene family with the observed sizes among taxa, 
CAFE uses a Monte Carlo re-sampling procedure. This option specifies the number of samples CAFE should use to calculate 
p-values. The tradeoff is between precision and computation time; in most cases 1000 samples should provide reasonable 
balance. Default: 1000.

\b -filter: The birth-death model of CAFE assumes at least one gene in the root of the species tree. This assumption may 
not be valid for families that were created after the most recent common ancestor of all species. The filter option 
filters out the families that are inferred (by parsimony) to have no genes in the root node of the species tree.
*/

#include <sstream>
#include <fstream>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <map>
#include <locale>
#include <set>

#include "gene_family.h"
#include <tree.h>

extern "C" {
#include "cafe.h"
extern pCafeParam cafe_param;
void __cafe_famliy_check_the_pattern(pCafeFamily pcf);
}

using namespace std;

class my_ctype : public std::ctype<char> {
  mask my_table[table_size + 1];
public:
  my_ctype(size_t refs = 0)
    : std::ctype<char>(my_table, false, refs)
  {
    // copy the original character classifaction table.
    std::copy(classic_table(),
      classic_table() + table_size,
      my_table);
    // and change ',' to be classed as white space.
    my_table[','] = (mask)space;
  }
};

vector<string> tokenize(string s, int flags)
{
  istringstream iss(s);
  if (flags & COMMA_AS_WHITESPACE)
  {
    std::locale x(std::locale::classic(), new my_ctype);
    iss.imbue(x);
  }
  vector<string> result;

  while (iss.good()) {
    string tmp;
    iss >> tmp;
    if (tmp.size() > 0)
      result.push_back(tmp);
  }

  return result;
}

vector<string> string_split(string line, char delimiter)
{
  vector<string> result;
  stringstream strstr(line);
  string word = "";
  while (getline(strstr, word, delimiter)) 
    result.push_back(word);
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

int to_int(std::string s)
{
	return stoi(s);
}

std::istream& operator>>(std::istream& ist, gene_family& fam)
{
	std::string str;
	std::getline(ist, str);
	std::vector<std::string> values = string_split(str, '\t');
	fam.desc = values[0];
	fam.id = values[1];
	fam.values.resize(values.size() - 2);
	std::transform(values.begin() + 2, values.end(), fam.values.begin(), to_int);
	return ist;
}


pCafeFamily load_gene_families(std::istream& ist, char separator, int max_size)
{
    char buf[STRING_BUF_SIZE];
    if (!ist)
    {
        return NULL;
    }

    ist.getline(buf, STRING_BUF_SIZE);
    string_pchar_chomp(buf);

    vector<string> species_list = string_split(buf, separator);
    if (species_list.size() < 2)
    {
        vector<string> t;
        species_list.swap(t); // guarantees that memory won't leak
        throw runtime_error("Failed to identify species for gene families");
    }
    species_list.erase(species_list.begin(), species_list.begin() + 2); // first two items are description and ID - delete them
    pCafeFamily pcf = cafe_family_init(species_list);

    string s;
    for (int i = 0; std::getline(ist, s); i++)
    {
        if (!ist)
            break;

        if (s.empty())
            continue;

        std::replace(s.begin(), s.end(), separator, '\t');
        istringstream isst(s);
        gene_family gf;
        isst >> gf;

        if (max_size < 0 || *max_element(gf.values.begin(), gf.values.end()) <= max_size)
            cafe_family_add_item(pcf, gf);
    }

    __cafe_famliy_check_the_pattern(pcf);

    return pcf;
}

/// Data array is assumes to contain a description, an identifier, and a set of integers
/// giving the family size in species order
void cafe_family_add_item(pCafeFamily pcf, const gene_family& gf)
{
  pCafeFamilyItem pitem = (pCafeFamilyItem)memory_new(1, sizeof(CafeFamilyItem));
  pitem->count = (int*)calloc(pcf->num_species, sizeof(int));
  pitem->maxlh = -1;   // Maximum likelihood index
  pitem->ref = -1;
  pitem->lambda = NULL;
  pitem->mu = NULL;
  pitem->holder = 1;

  if (gf.values.size() != size_t(pcf->num_species))
  {
    std::cerr << "Inconsistency in column count: expected " << pcf->num_species + 2 << ", but found " << gf.values.size() + 2;
  }

  pitem->desc = new char[gf.desc.size()+1];
  strcpy(pitem->desc, gf.desc.c_str());
  pitem->id = new char[gf.id.size() + 1];
  strcpy(pitem->id, gf.id.c_str());

  copy(gf.values.begin(), gf.values.end(), pitem->count);

  pcf->max_size = max(pcf->max_size, *std::max_element(pitem->count, pitem->count + pcf->num_species));
  arraylist_add(pcf->flist, pitem);
}

void cafe_family_item_free(pCafeFamilyItem pitem)
{
  delete [] pitem->id;
  pitem->id = NULL;
  memory_free(pitem->count);
  pitem->count = NULL;
  delete[] pitem->desc;
  pitem->desc = NULL;
  if (pitem->holder)
  {
    if (pitem->lambda) { memory_free(pitem->lambda); pitem->lambda = NULL; }
    if (pitem->mu) { memory_free(pitem->mu); pitem->mu = NULL; }
  }
  memory_free(pitem);
  pitem = NULL;
}


void cafe_family_filter(pCafeParam param)
{
  int i, n;
  pCafeFamily pcf = param->pfamily;
  pCafeTree pcafe = param->pcafe;
  pArrayList fflist = arraylist_new(11000);
  int max = 0;
  for (i = 0; i < pcf->flist->size; i++)
  {
    pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[i];
    cafe_family_set_size(pcf, pitem, pcafe);
    pArrayList nlist = pcafe->super.nlist;
    for (n = 0; n < nlist->size; n += 2)
    {
      pCafeNode pnode = (pCafeNode)nlist->array[n];
      if (pnode->familysize > 0)
      {
        pTreeNode p = ((pTreeNode)pnode);
        do
        {
          p->reg = 1;
          p = p->parent;
        } while (p && p->reg == 0);
      }
    }
    pTreeNode root = ((pTree)pcafe)->root;
    pTreeNode child[2] = { (pTreeNode)tree_get_child(root,0), (pTreeNode)tree_get_child(root,1) };

    root->reg = 0;
    if (child[0]->reg && child[1]->reg)
    {
      root->reg = 1;
    }

    // (1 and 1)
    else if (child[0]->reg || child[1]->reg)
    {
      root->reg = 0;
    }

    if (root->reg == 1)
    {
      arraylist_add(fflist, pitem);
      for (n = 0; n < pcf->num_species; n++)
      {
        if (pitem->count[n] > max)
        {
          max = pitem->count[n];
        }
      }
    }
    else
    {
      cafe_family_item_free(pitem);
    }
    tree_clear_reg((pTree)pcafe);
  }
  if (pcf->flist->size != fflist->size)
  {
    cafe_log(param, "The Number of families : %d ==> %d\n", pcf->flist->size, fflist->size);
    arraylist_free(pcf->flist, NULL);
    pcf->flist = fflist;
    for (i = 0; i < fflist->size; i++)
    {
      pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[i];
      pitem->ref = -1;
    }
    __cafe_famliy_check_the_pattern(pcf);
  }
  else
  {
    arraylist_free(fflist, NULL);
  }

  if (pcf->max_size != max)
  {
    pcf->max_size = max;
    init_family_size(&param->family_size, max);
    cafe_tree_set_parameters(pcafe, &param->family_size, 0);
  }
}

int log_cluster_membership(pCafeFamily pcf, int k_value, double **p_z_membership, std::ostream& log)
{
  if (p_z_membership == NULL) {
    fprintf(stderr, "family membership not found.\n");
    fprintf(stderr, "run lambdakmean or lambdamukmean to find cluster membership for each family.\n");
  }
  log << "The Number of families : " << pcf->flist->size << "\n";
  for (int i = 0; i < pcf->flist->size; i++)
  {
    pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[i];
    log << "family " << pitem->id << ": ";
    for (int k = 0; k < k_value; k++) {
      log << " " << p_z_membership[i][k];
    }
    log << endl;
  }
  return 0;
}

void cafe_family_reset_maxlh(pCafeFamily pcf)
{
  int i;
  pArrayList flist = pcf->flist;
  for (i = 0; i < flist->size; i++)
  {
    pCafeFamilyItem pitem = (pCafeFamilyItem)flist->array[i];
    pitem->maxlh = -1;
  }
}

int cafe_family_get_index(pCafeFamily pcf, const char* szid)
{
  int i;
  pArrayList flist = pcf->flist;
  for (i = 0; i < flist->size; i++)
  {
    pCafeFamilyItem pitem = (pCafeFamilyItem)flist->array[i];
    if (pitem->id && strcmp(pitem->id, szid) == 0) break;
  }
  return flist->size == i ? -1 : i;
}

string lower(char *cp)
{
    string n(cp);
    for (auto& c : n)
    {
        c = tolower(c);
    }
    return n;
}

/*! \brief Synchronize a family and tree together
*
*  Sets the value of index in pcf to the node index in the tree of the
*  leaf with the matching species ID.
*  \returns 0 on success, or -1 if there is a species in the tree with
*  no matching species in the family
*/void cafe_family_set_species_index(pCafeFamily pcf, pCafeTree pcafe)
{
    pTree ptree = (pTree)pcafe;
    map<string, int> leaf_names;

    for (int j = 0; j < ptree->nlist->size; j += 2)		// map leaf names to indices
    {
        pPhylogenyNode pnode = (pPhylogenyNode)ptree->nlist->array[j];
        leaf_names[lower(pnode->name)] = j;
    }

    set<string> all_species;
    for (int i = 0; i < pcf->num_species; i++)
    {
        string species = lower(pcf->species[i]);
        all_species.insert(species);
        if (pcf->species[i][0] == '-') {
            int nodeidx = atoi(&pcf->species[i][1]);
            pcf->index[i] = nodeidx;
        }
        else {
            if (leaf_names.find(species) == leaf_names.end())
                throw std::runtime_error("No species '" + species + "' was found in the tree");
            pcf->index[i] = leaf_names[species];
        }
    }

    for (auto& leaf : leaf_names)
    {
        if (all_species.find(leaf.first) == all_species.end())
            throw std::runtime_error("No species '" + leaf.first + "' was found in the family list");
    }
}
