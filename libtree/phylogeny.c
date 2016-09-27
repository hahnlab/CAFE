#include "tree.h"

#include<io.h>
#include<string.h>
#include<stdlib.h>
#include<sys/types.h>
#include<utils.h>
#include<utils_string.h>
#include<unistd.h>

void phylogeny_clear_node(pPhylogenyNode pnode)
{
	memset( pnode, 0, sizeof(PhylogenyNode) );
	pnode->branchlength = -1;
	pnode->bootstrap = -1;
	pnode->duplicated = PHYLOGENY_SPECIATION_NOT_DECIDED;
	pnode->taxaid = -1;
}

pTreeNode phylogeny_new_empty_node(pTree ptree)
{
	pPhylogenyNode pnode = (pPhylogenyNode) memory_new( 1, sizeof(PhylogenyNode) ); 
	phylogeny_clear_node(pnode);
	return (pTreeNode)pnode;
}

void phylogeny_node_copy(pTreeNode psrc, pTreeNode pdest )
{
	pPhylogenyNode ppsrc = (pPhylogenyNode)psrc;
	pPhylogenyNode ppdest = (pPhylogenyNode)pdest;

	ppdest->name = ppsrc->name;
	ppdest->gene = ppsrc->gene;
    ppdest->branchlength = ppsrc->branchlength;
	ppdest->species = ppsrc->species;
	ppdest->duplicated = ppsrc->duplicated;
	ppdest->taxaid = ppsrc->taxaid;
	ppdest->bootstrap = ppsrc->bootstrap;
}

typedef struct 
{
	char* name;
	char* branchlength;
	char* species;
	char* bootstrap;
	char* duplicated;
	char* gene;
	char* taxaid;
}PnodeParam;

typedef PnodeParam* pPnodeParam;


size_t __phylogeny_string_cmp(const void* arg1, const void* arg2)
{
	return strcasecmp((char*)arg1,(char*)arg2);
}

int phylogeny_set_node(pTree ptree, pTreeNode ptnode, pPnodeParam param )
{
	pPhylogenyNode pnode = (pPhylogenyNode)ptnode;
	if ( param->name )
	{
		pnode->name = (char*) memory_new(strlen(param->name)+1,sizeof(char));
		strcpy( pnode->name, param->name );
	}
	pnode->branchlength = -1;	
	if ( param->branchlength )
	{
		sscanf( param->branchlength, "%lf", &pnode->branchlength );
/*		if ( pnode->branchlength < 0 )
		{
			fprintf(stderr, "The branchlength of %s must be greater than 0\n", param->name);
			return -1;
		}
*/
	}
	if ( param->bootstrap )
	{
		pnode->bootstrap = atoi(param->bootstrap);
	}
	if ( param->duplicated )
	{
		pnode->duplicated = param->duplicated[0] == 'Y' ? 
				PHYLOGENY_SPECIATION_DUPLICATED : PHYLOGENY_SPECIATION_SPECIES;
	}
	else
	{
		pnode->duplicated = PHYLOGENY_SPECIATION_NOT_DECIDED;
	}

	if (param->taxaid )
	{
		sscanf( param->taxaid, "%d", &pnode->taxaid );
	}
	if ( param->gene )
	{
		pnode->gene = (char*) memory_new(strlen(param->gene)+1,sizeof(char));
		strcpy( pnode->gene, param->gene );
	}
	if (param->species)
	{
		char* sp = vector_get_by_cmp( (pVector)ptree->data, param->species, __phylogeny_string_cmp );
		if ( sp == NULL )
		{
			sp = (char*) memory_new(strlen(param->species)+1, sizeof(char));
			strcpy(sp,param->species);
			vector_add((pVector)ptree->data,sp);
		}		
		pnode->species = sp;
	}
	return 0;
}

void __phylogeny_free_node(pTree ptree, pTreeNode ptnode, va_list ap1)
{
	pPhylogenyNode pnode = (pPhylogenyNode)ptnode;
	if ( *ptree->count == 0  )
	{
		if ( pnode->name ) memory_free(pnode->name);
		pnode->name = NULL;
		if ( pnode->gene ) memory_free(pnode->gene);
		pnode->gene = NULL;
	}
	if ( ptnode->children ) vector_free( ptnode->children, NULL );
	if ( ptree->nlist ) 
	{
		ptree->nlist->array[ptnode->id] = NULL;
	}
	memory_free(ptnode);
	ptnode = NULL;
}

void phylogeny_free(pTree ptree)
{
	if ( *ptree->count == 0 && ptree->data ) vector_free( ((pVector)ptree->data), free );
	if ( ptree->nlist )
	{
		int i;
		for ( i = 0 ; i < ptree->nlist->size ; i++ )
		{
			__phylogeny_free_node(ptree, (pTreeNode)ptree->nlist->array[i], NULL );
		}
	}
	else
	{
		tree_traveral_prefix(ptree, __phylogeny_free_node);
	}
	tree_free(ptree);
}


pTree phylogeny_load_from_file(char* fname,
							   tree_func_new new_tree_func, 
							   tree_func_node_new new_tree_node_func, 
		                       phylogeny_func_parse_node parsefunc )
{
	int i;
	char* sztree = file_read_all(fname);	
	for( i = 0 ; i < sztree[i]; i++ )
	{
		if ( sztree[i] == '\n' ) sztree[i] = ' ';
	}
	pTree ptree = phylogeny_load_from_string(sztree, new_tree_func, 
							   new_tree_node_func, parsefunc ); 
	memory_free(sztree);
	sztree = NULL;
	return ptree;
}

char* phylogeny_interpret_node(pTree ptree, pTreeNode ptnode, char* sztree)
{
	register char c;
	int  state = 0;
	int  i = 0;
	PnodeParam param;
	memset( &param, 0, sizeof(PnodeParam) );
	while( *sztree )
	{ 	
		c = *sztree++;
		if ( c == '(' || c == ')' || c == ',' )
		{
			sztree--;
			break;
		}
		else
		{
			if ( state == 0 && c == ':' )
			{
				if ( i > 0 )
				{
					param.name = (char*)(sztree-i-1);
					*(char*)(sztree-1) = '\0';
				}
				i = 0;
				state = 1;
			}
			else if ( c == '[' )
			{
				if ( state == 0 && i >  0 )
				{
					param.name = (char*)(sztree-i-1);
					*(char*)(sztree-1) = '\0';
				}
				state = 2;

				if ( strncmp(sztree,"&&NHX",5) == 0 )
				{
					if ( i )
					{
						param.branchlength = (char*)(sztree-i-1);
						*(char*)(sztree-1) = '\0';
						i = 0;
					}
					char *next = index(sztree,']');	
					if ( next == NULL ) 
					{
						print_error(__FILE__,(char*)__FUNCTION__,__LINE__, "keep format : [&&NHX ... ]");		
					}
                    else {
                        *next++ = '\0';
                    }
					char* ptr = sztree + 5;
					char* tmp;
					while( *++ptr )
					{
						char* eq = index(ptr,'=');
						*eq++ = '\0';
						if ( (tmp= index(eq,':')))
						{
							*tmp = '\0';
						}
						if ( strcmp(ptr,"S") == 0 ) param.species = eq;
						else if ( strcmp(ptr,"B") == 0 ) param.bootstrap = eq;
						else if ( strcmp(ptr,"D") == 0 ) param.duplicated = eq;
						else if ( strcmp(ptr,"T") == 0 ) param.taxaid = eq;
						else if ( strcmp(ptr,"G") == 0 ) param.gene = eq;
						if ( tmp == NULL ) break;
						ptr = tmp;
					}
					sztree = next;
				}
			}
			else
			{
				i++;
			}
		}
	}
	char old = '\0';
	if ( i > 0 )
	{
		if ( state == 0 )
		{
			old = *sztree;
			param.name = (char*)(sztree-i);
			*sztree = '\0';
		}
		else if ( state == 1 )
		{
			old = *sztree;
			param.branchlength = (char*)(sztree-i);
			*sztree = '\0';
		}
	}

	if ( param.name )
	{
		size_t len = strlen(param.name);
		if ( len == 0 ) param.name = NULL;
		else if ( param.name[len-1] == '\n' && len == 1 ) param.name = NULL;
	}
	if ( param.branchlength )
	{
		size_t len = strlen(param.branchlength);
		if ( len == 0 ) param.branchlength = NULL;
		else if ( param.branchlength[len-1] == '\n' && len == 1 ) param.branchlength = NULL;
	}

	phylogeny_set_node(ptree, ptnode, &param);

	if ( old )
	{
		*sztree = old;
	}
	if ( *sztree == ',' )
	{
		sztree++;
	}
	return sztree;
}

int phylogeny_check_tree_string(char* sztree)
{
	size_t len = strlen(sztree);
	int i,grp_cnt = 0; 
	uintptr_t leaf_cnt = -1;
	int start = 0;
	int err = 0;
	pStack pstack = stack_new();
	for ( i = 0 ; i < len ; i++ )
	{
		if ( sztree[i] == ' ' ) continue;
		if ( sztree[i] == '(' ) 
		{
			start = i;
			if ( leaf_cnt != -1 ) stack_push( pstack, (void*)leaf_cnt );
			leaf_cnt = 0;
			grp_cnt++;
		}
		else if ( sztree[i] == ')') 
		{
/*			if ( leaf_cnt != 1 )
			{
				char buf[STRING_STEP_SIZE];
				strncpy( buf, &sztree[start], i - start + 1);
				buf[i-start+1] = '\0';
				fprintf(stderr,"Each node must contain two leaves. Check %s\n", buf );
				err = 2;
			}*/
			leaf_cnt = (uintptr_t)stack_pop( pstack );
			grp_cnt--;
		}
		else if ( sztree[i] == ',' )
		{
			leaf_cnt++;
		}
	}
	if ( grp_cnt != 0 )
	{
		fprintf(stderr, "Tree error (Unbalanced parentheses): %s\n", sztree);
		err = 1;
	}	
	stack_free(pstack);
	return err;
}

pTree phylogeny_load_from_string(char* sztree, 
								 tree_func_new new_tree_func, 
								 tree_func_node_new new_tree_node_func, 
		                         phylogeny_func_parse_node parsefunc, ... )
{
	va_list ap; 
	va_start(ap,parsefunc);
	pTree ptree = new_tree_func(new_tree_node_func, ap);
	va_end(ap);
	ptree->data = vector_new();
	pTreeNode cur = ptree->root;
	int bfirst = 1;
	char c;
	if ( phylogeny_check_tree_string(sztree) )
	{
		return NULL;
	}

	string_pchar_chomp(sztree);
	size_t tlen = strlen(sztree);
	if ( sztree[tlen-1] == ';' )
	{
		sztree[tlen-1] = '\0';
	}

	while( (c = *sztree++) )
	{
		if ( c == ' ' || c == '\t' || c == '\n' || c == '\r' ) continue;
		if ( c == ';' ) break;

		if ( c == '(' )
		{
			if ( bfirst )
			{
				bfirst = 0;
				continue;
			}
			pTreeNode ptnode = new_tree_node_func(ptree);
			if ( cur->children == NULL )
			{
				cur->children = vector_new();
			}
			vector_add(cur->children, ptnode);
			ptnode->parent = cur;
			cur = ptnode;
		}
		else if ( c == ')' )
		{
			if ( !*sztree ) break;
			sztree = phylogeny_interpret_node(ptree,cur,sztree);
			if ( parsefunc) parsefunc(ptree,cur);
			cur = cur->parent;
		}
		else
		{
			sztree--;
			pTreeNode ptnode = new_tree_node_func(ptree);
			ptnode->parent = cur;
			sztree = phylogeny_interpret_node(ptree,ptnode,sztree);
			if ( parsefunc ) parsefunc(ptree,ptnode);
			if ( cur->children == NULL )
			{
				cur->children = vector_new();
			}
			vector_add(cur->children, ptnode);
		}
	}
	return ptree;
}

pTree phylogeny_root_tree(pTree ptree, 
							tree_func_node_new new_tree_node_func, 
							phylogeny_func_parse_node parsefunc, ... )
{
	pTreeNode cur = ptree->root;	
	pTreeNode dummynode = new_tree_node_func(ptree);
	if ( dummynode->children == NULL )
	{
		dummynode->children = vector_new();
	}
	while (cur->children->size > 1) {
		pTreeNode child = tree_get_child(cur, 0);
		vector_dereference_by_data(cur->children, child);	
		tree_add_child_to_node(dummynode, child);
	}
	pTreeNode child = tree_get_child(cur, 0);
	vector_dereference_by_data(cur->children, child);	
	tree_add_child_to_node(cur, dummynode);
	tree_add_child_to_node(cur, child);
	return ptree;
}

void __phylogeny_string(pString pstr, pPhylogenyNode pnode, phylogeny_func_name_modify fmod )
{
	char buf[STRING_BUF_SIZE];
	int reg = pnode->super.reg;
	if ( pnode->name || fmod )
	{
		if ( fmod )
		{
			fmod(pstr, pnode);
		}
		else
		{
			string_add(pstr,pnode->name);				
		}
	}
	if ( !(reg&PS_SKIP_BL) && pnode->branchlength >= 0 )
	{
		sprintf(buf,":%g", pnode->branchlength);
		string_add(pstr,buf);
	}
	if ( !(reg & PS_NHX ) ) return;
	if ( pnode->species || pnode->bootstrap > 0 || 
		 pnode->duplicated != PHYLOGENY_SPECIATION_NOT_DECIDED ||
		 pnode->taxaid != -1 )
	{
		string_add(pstr,"[&&NHX");
		if ( pnode->species )
		{
			sprintf(buf,":S=%s", pnode->species);
			string_add(pstr,buf);
		}
		if ( pnode->bootstrap >= 0 )
		{
			sprintf(buf,":B=%d", pnode->bootstrap);
			string_add(pstr,buf);
		}
		if ( pnode->duplicated != PHYLOGENY_SPECIATION_NOT_DECIDED )
		{
			sprintf(buf,":D=%s", pnode->duplicated == PHYLOGENY_SPECIATION_DUPLICATED  ? "Y" : "N" );
			string_add(pstr,buf);
		}
		/*
		if ( pnode->gene )
		{
			sprintf(buf,":G=%s", pnode->gene );
			string_add(pstr,buf);
		}
		*/
		if ( pnode->taxaid != -1 )
		{
			sprintf(buf,":T=%d", pnode->taxaid );
			string_add(pstr,buf);
		}
		string_add(pstr,"]");
	}
}

pString phylogeny_string_newick(pTree ptree, phylogeny_func_name_modify fmod, int opts)
{
	int i;
	pArrayList nlist = ptree->nlist;
	for ( i = 0 ; i < nlist->size ; i++ )
	{
		pTreeNode pnode = (pTreeNode)nlist->array[i];	
		pnode->reg = opts;
	}
	return phylogeny_string(ptree, fmod );
}

pString phylogeny_string_nhx(pTree ptree, phylogeny_func_name_modify fmod, int opts)
{
	int i;
	pArrayList nlist = ptree->nlist;
	for ( i = 0 ; i < nlist->size ; i++ )
	{
		pTreeNode pnode = (pTreeNode)nlist->array[i];	
		pnode->reg = (PS_NHX & opts) | 2;
	}
	return phylogeny_string(ptree, fmod );
}


pString phylogeny_string(pTree ptree, phylogeny_func_name_modify fmod)
{
	pString pstr = string_new();
	pStack pstack = stack_new();
	stack_push(pstack,ptree->root);	
	while( stack_is_empty(pstack) )
	{
		pTreeNode ptnode = (pTreeNode)pstack->head->data;
		pPhylogenyNode pnode = (pPhylogenyNode)ptnode;

		if ( ptnode->children )
		{
			pLinkedListItem child = ptnode->children->tail;
			if ( child && (((pTreeNode)child->data)->reg & 1) )
			{
				string_add(pstr,")");
				stack_pop(pstack);
				__phylogeny_string(pstr, pnode, fmod);
				ptnode->reg |= 1;
			}
			else
			{
				if ( ptnode != ptree->root && 
				     ptnode->parent->children->size == 2 &&
			    	((pTreeNode)ptnode->parent->children->head->data)->reg & 1 )
				{
					string_add(pstr,",");
				}
				if ( ptnode->children->size == 2 )
				{
					string_add(pstr,"(");
				}
				while( child )
				{
					stack_push(pstack,child->data);	
					child = child->prev;			
				}
			}
		}
		else
		{
			stack_pop(pstack);
			if ( ptnode != ptree->root && 
			    ptnode->parent->children->size == 2 &&
			    ((pTreeNode)ptnode->parent->children->head->data)->reg & 1 )
			{
				string_add(pstr,",");
			}
			__phylogeny_string(pstr, pnode, fmod);
			ptnode->reg |= 1;
		}
	}
	stack_free(pstack);
	tree_clear_reg(ptree);
	return pstr;
}


pTree phylogeny_new(char* sztree, phylogeny_func_parse_node parsefunc )
{
	return	phylogeny_load_from_string(sztree, tree_new, phylogeny_new_empty_node, parsefunc );
}

pTree phylogeny_copy(pTree psrc)
{
	return tree_copy(psrc, phylogeny_new_empty_node, phylogeny_node_copy );
}

pTreeNode phylogeny_get_sibling(pTree ptree, pTreeNode ptnode )
{
	if ( tree_is_root(ptree, ptnode) ) return NULL;
	pTreeNode parent = ptnode->parent;					
	pTreeNode sibling = NULL;
	vector_rewind(parent->children);
	while( (sibling= (pTreeNode)vector_next(parent->children)) )
	{
		if ( sibling != ptnode ) break;
	}
	return sibling;
}

pTree phylogeny_split_tree(pTree ptree, int idx, tree_func_node freenode ) 
{
	pTreeNode pnode = (pTreeNode)ptree->nlist->array[idx];
	if ( tree_is_root(ptree,(pTreeNode)pnode) ) return NULL;			
	pTree psubtree = (pTree)memory_new(1,ptree->size);
//	memcpy(psubtree,ptree, sizeof(ptree->size));
	psubtree->count = ptree->count;
	(*psubtree->count)++;
	psubtree->size = ptree->size;
	psubtree->nlist = NULL;
	psubtree->postfix = NULL;
	psubtree->prefix = NULL;

	psubtree->root = pnode;
	pTreeNode sibling = phylogeny_get_sibling(ptree,pnode);
	if ( tree_is_root(ptree,pnode->parent) ) 
	{
		freenode(ptree, ptree->root,NULL);	
		ptree->root = sibling;	
		ptree->root->parent = NULL;	
	}
	else
	{
		((pPhylogenyNode)sibling)->branchlength += ((pPhylogenyNode)sibling->parent)->branchlength;
		pTreeNode parent = pnode->parent;
		pTreeNode grand = parent->parent;
		if ( grand->children->head->data == parent )
		{
			grand->children->head->data = sibling;
		}
		else
		{
			grand->children->tail->data = sibling;
		}
		sibling->parent = grand;
		freenode(ptree,parent,NULL);	
	}
	tree_build_node_list(ptree);
	tree_build_node_list(psubtree);
	psubtree->root->parent = NULL;
	((pPhylogenyNode)ptree->root)->branchlength = -1;
	((pPhylogenyNode)psubtree->root)->branchlength = -1;
	return psubtree;
}

void __phylogeny_delete_nodes(pTree ptree, pTreeNode pnode, tree_func_node freenode ) 
{
	pTreeNode sibling = phylogeny_get_sibling(ptree,pnode);
	if ( tree_is_root(ptree,pnode->parent) ) 
	{
		freenode(ptree, ptree->root,NULL);	
		ptree->root = sibling;	
		ptree->root->parent = NULL;	
	}
	else
	{
		((pPhylogenyNode)sibling)->branchlength += ((pPhylogenyNode)sibling->parent)->branchlength;
		pTreeNode parent = pnode->parent;
		pTreeNode grand = parent->parent;
		if ( grand->children->head->data == parent )
		{
			grand->children->head->data = sibling;
		}
		else
		{
			grand->children->tail->data = sibling;
		}
		sibling->parent = grand;
		freenode(ptree,parent,NULL);	
	}

	Tree tmp;
	pnode->parent = NULL;
	tmp.nlist = ptree->nlist;
	tmp.root = pnode;
	tmp.count = ptree->count;
	tree_traveral_prefix(&tmp, freenode );
}

int phylogeny_delete_nodes_by_index(pTree ptree, int idx, tree_func_node freenode ) 
{
	pTreeNode pnode = (pTreeNode)ptree->nlist->array[idx];
	if ( tree_is_root( ptree, pnode) )
	{
		phylogeny_free(ptree); 
		return -1;
	}
	__phylogeny_delete_nodes(ptree, pnode, freenode );
	tree_build_node_list(ptree);
	((pPhylogenyNode)ptree->root)->branchlength = -1;
	return 0;
}

int phylogeny_delete_nodes_by_reg(pTree ptree, tree_func_node freenode ) 
/*
 * reg == 0 : delete
 * reg == 1 : remain
 */
{
	if ( ptree->root->reg == 0 )
	{
		phylogeny_free(ptree); 
		return -1;
	}
	int i, bupdate = 0;
	pArrayList nlist = ptree->nlist;

	for ( i = 0 ; i < nlist->size ; i++ )
	{
		pTreeNode pnode = (pTreeNode)nlist->array[i];
		if ( pnode == NULL ) continue;
		if ( ptree->root == pnode ) continue;
		// Find the upper most node to be deleted
		if ( pnode->reg == 0  ) 
		{
			pTreeNode p;
			for ( p = pnode->parent ; p ; p = p->parent )
			{
				if ( p->reg == 1 ) break;
				pnode = p;
			}
			if ( p ) 
			{
				__phylogeny_delete_nodes(ptree,pnode,freenode);
				bupdate = 1;
			}
		}
	}
	for ( i = 0 ;  i < nlist->size ; i++ ) 
	{
		pTreeNode p = nlist->array[i];
		if ( p ) p->reg = 0;
	}
	if ( bupdate )
	{
		((pPhylogenyNode)ptree->root)->branchlength = -1;
		tree_build_node_list(ptree);
		/*
		pPhylogenyNode root = (pPhylogenyNode)ptree->root;
		if ( root->duplicated != PHYLOGENY_SPECIATION_NOT_DECIDED )
		{
			root->duplicated = PHYLOGENY_SPECIATION_SPECIES;
		}
		*/
	}
	tree_clear_reg((pTree)ptree);
	return 0;
}

int phylogeny_delete_nodes_by_func(pTree ptree, tree_func_node filter ,tree_func_node freenode ) 
{
	int i;
	pArrayList nlist = ptree->nlist;
	for ( i = 0 ; i < nlist->size ; i+=2 )
	{
		pTreeNode pnode = (pTreeNode)nlist->array[i];
		filter(ptree,pnode,NULL);
		if ( pnode->reg )
		{
			while( pnode && pnode->reg == 0 )
			{
				pnode->reg = 1;			
				pnode = pnode->parent;
			}
		}
	}
	return phylogeny_delete_nodes_by_reg(ptree, freenode );
}

void phylogeny_increase_branchlength(pTree ptree, int idx)
{	
	pPhylogenyNode pnode = (pPhylogenyNode)ptree->nlist->array[idx];
	pnode->branchlength++;
}

void phylogeny_add_branchlength(pTree ptree, int idx, int add)
{
	pPhylogenyNode pnode = (pPhylogenyNode)ptree->nlist->array[idx];
	pnode->branchlength+=add;
}

pString phylogeny_to_mp(pTree ptree, pMetapostConfig pmc, ... )
{
	va_list ap;
	va_start(ap,pmc);

	pString pstr = string_new();
	pArrayList nlist = ptree->nlist;
	double* gloc = (double*) memory_new( nlist->size, sizeof(double));

	double gunit = 0;
	switch( pmc->dir )
	{
		case MP_DIR_VERTICAL:
			gunit = (double)pmc->width/(double)((nlist->size+1)/2 - 1);
			break;
		default:
			gunit = (double)pmc->height/(double)((nlist->size+1)/2 - 1);
			break;
	}

	double nunit = 0;
	int i, depth;
	for ( i = 0 ; i < nlist->size; i+=2 )
	{
		gloc[i] = gunit*(i>>1);
		pTreeNode pnode = (pTreeNode)nlist->array[i];
		pnode->reg = 0;
		pnode = pnode->parent;
		depth = 0;
		while( pnode )
		{
			depth++;
			if( pnode->reg > depth ) break;
			pnode->reg = depth;
			pnode = pnode->parent;
		}
	}
	if ( pmc->dir != MP_DIR_VERTICAL )
	{
		for ( i = 0 ; i < nlist->size; i+=2 )
		{
			gloc[i] = pmc->height - gloc[i];
		}
	}

	depth = ptree->root->reg;
	switch( pmc->dir )
	{
		case MP_DIR_VERTICAL:
			nunit = (double)pmc->height/depth;
			break;
		default:
			nunit = (double)pmc->width/depth;
			break;
	}

	if ( pmc->shape & MP_SHAPE_MOST_CENTER )	
	{
		for( i = 1 ; i < nlist->size; i+=2 )
		{
			pTreeNode pnode = (pTreeNode)nlist->array[i];		
			pTreeNode left, right;
			left = right = NULL;
	
			while( pnode->children )
			{
				left = (pTreeNode)pnode->children->head->data;
				pnode = left;
			}
	
			pnode = (pTreeNode)nlist->array[i];		
			while( pnode->children )
			{
				right = (pTreeNode)pnode->children->tail->data;
				pnode = right;
			}
			gloc[i] = (gloc[left->id] + gloc[right->id])/2.0;
		}
	}
	else if ( pmc->shape & MP_SHAPE_LEAST_CENTER )
	{
		for( i = 1 ; i < nlist->size; i+=2 )
		{
			gloc[i] = (gloc[i-1] + gloc[i+1])/2.0;
		}
	}

	string_fadd( pstr, "beginfig(%d)\n", pmc->id );
	switch( pmc->unit )
	{
		case MP_UNIT_MM: string_add( pstr, "u:=1mm;\n"); break;
		case MP_UNIT_CM: string_add( pstr, "u:=1cm;\n"); break;
		case MP_UNIT_IN: string_add( pstr, "u:=1in;\n"); break;
	}

// Point
	string_add( pstr, "pair p[];\n");
	string_add( pstr, "pair mid[];\n");
	for ( i = 0 ; i < nlist->size ; i++ )
	{
		pTreeNode pnode = (pTreeNode)nlist->array[i];
		switch( pmc->dir )
		{
			case MP_DIR_VERTICAL:
				string_fadd( pstr, "p%d = (%fu,%fu);\n", i, gloc[i], pnode->reg*nunit );

				break;
			default:
				string_fadd( pstr, "p%d = (%fu,%fu);\n", i, (depth-pnode->reg)*nunit, gloc[i] );
				break;
		}
	}
	string_add( pstr, "\n");
	tree_clear_reg(ptree);

// Draw line	
	string_add( pstr, "\% draw line\n");


	for ( i = 0; i < nlist->size ; i++ )
	{
		pTreeNode pnode = (pTreeNode)nlist->array[i];
		if ( tree_is_root( ptree, pnode )) 
		{
			if ( (pmc->shape & 3) != MP_SHAPE_TRI ) 
			{
				switch( pmc->dir )
				{
					case MP_DIR_VERTICAL:
						string_fadd( pstr, "draw p[%d] -- ( p[%d] + (0,%fu) )", i, i, nunit/7.0 ); 
						break;
					default:
						string_fadd( pstr, "draw p[%d] -- ( p[%d] + (-%fu,0) )", i, i, nunit/7.0 ); 
						break;
				}
			}
		}
		else
		{
			switch ( pmc->shape & 3 )
			{
				case MP_SHAPE_RECT:
					switch( pmc->dir )
					{
						case MP_DIR_VERTICAL:
				//			string_add( pstr, "\tdraw a -- ( xpart a, ypart b ) -- b ");
							string_fadd( pstr, "draw p[%d] -- ( xpart p[%d], ypart p[%d] ) -- p[%d] ", 
										i, i, pnode->parent->id, pnode->parent->id);
							break;
						default:
				//			string_add( pstr, "\tdraw a -- ( xpart b, ypart a ) -- b ");
							string_fadd( pstr, "draw p[%d] -- ( xpart p[%d], ypart p[%d] ) -- p[%d] ", 
									    i, pnode->parent->id, i, pnode->parent->id);
							break;
					}
					break;
				case MP_SHAPE_TRI:
					string_fadd(pstr, "draw p[%d] -- p[%d] ", i, pnode->parent->id );
					break;
			}
		}
		if ( pmc->fmod ) pmc->fmod(pstr,pnode, pmc, ap);
		else string_add(pstr,";\n");
	}

// Add name
	string_add( pstr, "\n\% add name\n");
	for ( i = 0 ; i < nlist->size ; i+=2)
	{
		pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
		switch( pmc->dir )
		{
			case MP_DIR_VERTICAL:
			//	string_fadd( pstr, "dotlabel.bot(btex \\small{%s} etex rotated 90, p[%d]);\n", pnode->name, i);
				string_fadd( pstr, "labeloffset := 10; dotlabel.bot(btex %s etex, p[%d]);\n", pnode->name, i);
				break;
			default:
				string_fadd( pstr, "labeloffset := 10; dotlabel.rt(btex %s etex, p[%d]);\n", pnode->name, i);
				break;
		}
	}

// Duplicate
	string_add( pstr, "\n\% Mark duplication\n");
	for ( i = 1 ; i < nlist->size ; i+=2)
	{
		pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
		if ( pnode->duplicated == PHYLOGENY_SPECIATION_DUPLICATED )
		{
			string_fadd( pstr, "draw p[%d] withpen pencircle scaled 4 withcolor red;\n", i );
		}
	}

	if ( pmc->fremark ) 
	{
		pmc->fremark(pstr,ptree,pmc,ap);
	}
	string_fadd( pstr, "endfig;\n");
	memory_free(gloc);
	gloc = NULL;
	va_end(ap);
	tree_clear_reg(ptree);
	return pstr;
}
