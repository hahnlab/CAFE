#ifndef __TREE_H__
#define __TREE_H__

#include<utils.h>
#include<stdarg.h>
#include<utils_string.h>

typedef struct tagTree* pTree;
typedef struct tagTreeNode* pTreeNode;

typedef struct tagTree
{
	pTreeNode 	 root;
	void*		 data;
	int			 size;
	int*		 count;
	pArrayList 	 nlist;
	pArrayList   postfix, prefix;
}Tree;

typedef struct tagTreeNode 
{
	pTreeNode 	parent;
	pVector 	children;
	int			reg;
	int 		id;
}TreeNode;

/****************************************************************************
 * Tree
****************************************************************************/

typedef pTreeNode (*tree_func_node_new)(pTree ptree);
typedef pTree (*tree_func_new)(tree_func_node_new nfunc, va_list ap);
typedef void (*tree_func_node)(pTree ptree, pTreeNode pnode, va_list ap);
typedef void (*tree_func_node_copy)(pTreeNode psrc, pTreeNode pdest );

extern pTree tree_new(tree_func_node_new nfunc, va_list ap);
extern void tree_new_fill(pTree ptree, tree_func_node_new nfunc ) ;
extern void tree_free(pTree ptree);
extern pTreeNode tree_new_empty_node();
extern int tree_is_leaf(pTreeNode pnode);
extern int tree_is_root(pTree ptree, pTreeNode pnode);
extern void tree_add_child_to_node(pTreeNode parent, pTreeNode child);
extern void tree_traveral_prefix(pTree ptree, tree_func_node func, ...);
extern void tree_traveral_postfix(pTree ptree, tree_func_node func, ...);
extern void tree_traveral_infix(pTree ptree, tree_func_node func, ... );
extern void tree_clear_reg(pTree ptree);
extern void* tree_get_child(pTreeNode ptnode, int idx);
extern pTree tree_copy(pTree psrc, tree_func_node_new new_node_func, tree_func_node_copy copy);
extern void tree_build_node_list(pTree ptree);

typedef enum
{
	PHYLOGENY_SPECIATION_NOT_DECIDED = 0,
	PHYLOGENY_SPECIATION_DUPLICATED,
	PHYLOGENY_SPECIATION_SPECIES 
}enumPhylogeny_speciation ;

typedef struct tagPhylogenyNode* pPhylogenyNode;
typedef struct tagPhylogenyNode
{
	TreeNode super;
	char*	name;	
	double  branchlength;
	char*	species;
	char*   gene;
	int 	duplicated;
	int		bootstrap;
	int 	taxaid;
}PhylogenyNode;

/****************************************************************************
 * Phylogeny
****************************************************************************/

typedef void (*phylogeny_func_parse_node)(pTree ptree, pTreeNode pnode);
typedef void (*phylogeny_func_name_modify)(pString str, pPhylogenyNode pnode);

extern pTree phylogeny_split_tree(pTree ptree, int idx, tree_func_node freenode );
extern pString phylogeny_string(pTree ptree, phylogeny_func_name_modify fmod);
extern void phylogeny_node_copy(pTreeNode psrc, pTreeNode pdest );
extern pTree phylogeny_load_from_string(char* sztree, 
								 tree_func_new new_tree_func, 
								 tree_func_node_new new_tree_node_func, 
		                         phylogeny_func_parse_node parsefunc, ... );
extern pTree phylogeny_root_tree(pTree ptree, 
							tree_func_node_new new_tree_node_func, 
							phylogeny_func_parse_node parsefunc, ... );
extern pTree phylogeny_load_from_file(char* sztree, 
								 tree_func_new new_tree_func, 
								 tree_func_node_new new_tree_node_func, 
		                         phylogeny_func_parse_node parsefunc );
extern void phylogeny_free(pTree ptree);

extern pTreeNode phylogeny_new_empty_node(pTree ptree);
extern void phylogeny_clear_node(pPhylogenyNode pnode);
extern pTree phylogeny_copy(pTree psrc);
extern pTree phylogeny_new(char* sztree, phylogeny_func_parse_node parsefunc );
extern int phylogeny_delete_nodes_by_index(pTree ptree, int idx, tree_func_node freenode );
extern int phylogeny_delete_nodes_by_reg(pTree ptree, tree_func_node freenode );
extern int phylogeny_delete_nodes_by_func(pTree ptree, tree_func_node filter ,tree_func_node freenode );
extern pTreeNode phylogeny_get_sibling(pTree ptree, pTreeNode ptnode );

extern pString phylogeny_string_nhx(pTree ptree, phylogeny_func_name_modify fmod, int opts);
extern pString phylogeny_string_newick(pTree ptree, phylogeny_func_name_modify fmod, int opts);


#define PS_SKIP_BL	0x1000
#define PS_NWICK	0x0000
#define PS_NHX 		0x0FFE


/****************************************************************************
 * Gene and Species Tree
****************************************************************************/
typedef struct
{
	PhylogenyNode super;
	pVector vsubspecies;
	int* 	mask;
	int 	gain;
	int 	numduplicated;
	int 	size;			// family size
} SpeciesTreeNode;

typedef SpeciesTreeNode* pSpeciesTreeNode;

typedef struct
{
	PhylogenyNode super;
	pVector vsubspecies;
	pSpeciesTreeNode psnode;
}GeneTreeNode;

typedef struct
{
	Tree super;
	pTree species_tree;
	int* species_tid;
}GeneTree;

typedef GeneTree* pGeneTree;

typedef GeneTreeNode* pGeneTreeNode;

extern void gene_tree_free(pGeneTree ptree);
extern void species_tree_free(pTree ptree);
extern pTree species_tree_new(char* sztree, phylogeny_func_parse_node fparse );
extern pGeneTree gene_tree_new(char* sztree, pTree pstree, phylogeny_func_parse_node fparse );
extern pString gene_tree_string(pGeneTree ptree, int mode);
extern pString species_tree_string(pTree ptree, int mode);
extern int gene_tree_map_species(pGeneTree pgtree);
extern int gene_tree_duplication_count(pGeneTree pgtree);
extern void gene_tree_duplication_count_report(pTree pstree, FILE* fout);
extern void species_tree_clear_gain(pTree pstree);

extern pString species_tree_metapost(pTree pstree, int id, char* title, double width, double height );
extern pString gene_tree_metapost(pGeneTree pstree, int id, char* title, double width, double height );

/****************************************************************************
 * MetaPhost Config
****************************************************************************/
typedef enum
{
	MP_UNIT_IN = 0,
	MP_UNIT_MM, 
	MP_UNIT_CM
}tagMPUnitList;

typedef enum
{
	MP_DIR_HORIZONTAL = 0,
	MP_DIR_VERTICAL 
}tagMPDirectionList;

typedef enum
{
	MP_SHAPE_RECT = 1,
	MP_SHAPE_TRI = 2,
	MP_SHAPE_MASK = 3,
	MP_SHAPE_MOST_CENTER = 4,
	MP_SHAPE_LEAST_CENTER = 8,
}tagMPShapeList;

typedef struct tagMetapostConfig* pMetapostConfig;
typedef double (*metapost_remark)(pString str, pTree ptree, pMetapostConfig  pmc, va_list ap);
typedef double (*metapost_annotation_func)(pString str, pTreeNode pnode, pMetapostConfig pmc, va_list ap);

typedef struct tagMetapostConfig
{
	int id;
	int unit;
	int dir;
	int shape;
	metapost_remark fremark;
	metapost_annotation_func fmod;
	double width, height;
}MetapostConfig;


pString phylogeny_to_mp(pTree ptree, pMetapostConfig pmc, ... );

#define TREE_DEBUG 1

#endif
