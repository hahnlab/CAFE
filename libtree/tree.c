#include "tree.h"
#include<fcntl.h>
#include<string.h>
#include<stdlib.h>

void tree_new_fill(pTree ptree, tree_func_node_new nfunc ) 
{
	ptree->size = sizeof(Tree);
	if ( nfunc ) ptree->root = nfunc(ptree);
	ptree->count = (int*)memory_new(1,sizeof(int));
	*ptree->count = 0;
	ptree->nlist = NULL;
	ptree->postfix = NULL;
	ptree->prefix = NULL;
}

pTree tree_new(tree_func_node_new nfunc, int size) 
{
	pTree ptree = (pTree)memory_new(1,sizeof(Tree));
	if ( nfunc == NULL ) 
	{
		nfunc = tree_new_empty_node;
	}
	tree_new_fill(ptree, nfunc);
	return ptree;
}

void tree_free(pTree ptree)
{
	if ( *ptree->count == 0 )
	{
		memory_free(ptree->count);
		ptree->count = NULL;
	}
	else
	{
		(*ptree->count)--;
	}
	if ( ptree->nlist )
	{
		arraylist_free(ptree->nlist,NULL);		
		arraylist_free(ptree->postfix,NULL);		
		arraylist_free(ptree->prefix,NULL);		
	}
	memory_free(ptree);
	ptree = NULL;
}

pTreeNode tree_new_empty_node(pTree ptree)
{
	pTreeNode ptnode = (pTreeNode)memory_new(1, sizeof(TreeNode));
	ptnode->reg = 0;
	ptnode->children = NULL;
	ptnode->parent = NULL;
	return ptnode;
}

int tree_is_leaf(pTreeNode pnode)
{
	return pnode->children == NULL;
}

int tree_is_root(pTree ptree, pTreeNode pnode)
{
	return ptree->root == pnode;
}

void tree_add_child_to_node(pTreeNode parent, pTreeNode child)
{
	child->parent = parent;			
	if ( parent->children == NULL )
	{
		parent->children = vector_new();
	}
	vector_add(parent->children, child);
}


void* tree_get_child(pTreeNode ptnode, int idx)
{
	if ( ptnode->children == NULL ) return NULL;			
	return vector_get(ptnode->children, idx);
}

int tree_get_child_index(pTree ptree, pTreeNode ptnode)
{
	if ( tree_is_root(ptree,ptnode) ) return -1;
	pTreeNode parent = ptnode->parent;
	pLinkedListItem cur = parent->children->head;
	int i;
	for( i = 0 ; cur ; i++, cur = cur->next )
	{
		if ( cur->data == ptnode ) break;
	}
	return i;
}

// Only for binary tree

void tree_traveral_prefix(pTree ptree, tree_func_node func, ...)
{
	va_list ap;
	va_start( ap, func );
	pStack pstack = stack_new();
	stack_push(pstack,ptree->root);	
	while( stack_has_items(pstack) )
	{
		pTreeNode pnode = (pTreeNode)stack_pop(pstack);
		// from right to left
		if ( pnode->children )
		{
			pLinkedListItem child = pnode->children->tail; 
			while( child )
			{
				stack_push(pstack,child->data);	
				child = child->prev;
			}
		}
		func(ptree, pnode, ap);
	}
	va_end(ap);
	stack_free(pstack);
}

void tree_traveral_postfix(pTree ptree, tree_func_node func, ...)
{
	va_list ap;
	va_start( ap, func );
	pStack pstack = stack_new();
	stack_push(pstack,ptree->root);	
	while( stack_has_items(pstack) )
	{
		pTreeNode pnode = (pTreeNode)pstack->head->data;
		int visit_count = -1;
		if ( pnode->children )
		{
			visit_count = 0;
			pLinkedListItem child = pnode->children->tail; 
			while( child )
			{
				if ( ((pTreeNode)child->data)->reg )
				{
					visit_count++;
				}
				else
				{
					stack_push(pstack,child->data);	
				}
				child = child->prev;			
			}
			if ( visit_count == pnode->children->size )
			{
				visit_count = -2;
			}
		}
		if ( visit_count < 0 )
		{
			stack_pop(pstack);
			func(ptree, pnode, ap);
			pnode->reg = 1;
		}
	}
	va_end(ap);
	stack_free(pstack);
	tree_clear_reg(ptree);
}

// Only for complete binary tree
void tree_traveral_infix(pTree ptree, tree_func_node func, ... )
{
	va_list ap;
	va_start( ap, func );
	pStack pstack = stack_new();
	stack_push(pstack,ptree->root);	
	while( stack_has_items(pstack) )
	{
		pTreeNode pnode = (pTreeNode)pstack->head->data;
		if ( pnode->children )
		{
			if ( pnode->children->size != 2 ) {
				fprintf(stderr, "ERROR(tree.c): running infix on nonbinary tree\n");
				return;
			} 
			pTreeNode left = (pTreeNode)pnode->children->head->data;
			if ( left->reg )
			{
				stack_pop(pstack);
				pnode->reg = 1;
				if (  pnode->children->size > 1  )
				{
					stack_push(pstack, pnode->children->tail->data);
				}
				func(ptree, pnode, ap);
			}
			else
			{
				stack_push(pstack, left);	
			}
		}
		else
		{
			stack_pop(pstack);
			func(ptree, pnode, ap);	// adds to nlist
			pnode->reg = 1;
		}
	}
	va_end(ap);
	stack_free(pstack);
	tree_clear_reg(ptree);
}


void __tree_clear_reg(pTree ptree, pTreeNode pnode, va_list ap1)
{
#ifdef __DEBUG_TREE__
	if ( tree_is_leaf(pnode) )
	{
		fprintf(stderr, ">> %s\n", ((pPhylogenyNode)pnode)->name );
	}
	else
	{
		if ( ptree->root == pnode )
		{
			fprintf(stderr, "Root\n");
		}
		else
		{
			fprintf(stderr, "Inter\n");
		}
	}
#endif
	pnode->reg = 0;
}

void tree_clear_reg(pTree ptree)
{
	if ( ptree->nlist )
	{
		int i;
		for ( i = 0 ; i < ptree->nlist->size ; i++ )
		{
			((pTreeNode)ptree->nlist->array[i])->reg = 0;
		}
	}
	else
	{
		tree_traveral_prefix(ptree, __tree_clear_reg );
	}
}

void __tree_copy(pTree psrc, pTreeNode pnode, va_list ap1)
{
  va_list ap;
  va_copy(ap, ap1);
	pTree pdest = va_arg(ap,pTree);
	tree_func_node_new new_node_func = va_arg(ap,tree_func_node_new);
	tree_func_node_copy copy = va_arg(ap,tree_func_node_copy);
	pTreeNode cur = (pTreeNode)stack_pop((pStack)pdest->data);
	copy(pnode,cur);
	if ( !tree_is_leaf(pnode) )
	{
		int i;
		cur->children = vector_new();
		for ( i = 0 ; i < pnode->children->size; i++ )
		{
			pTreeNode child = new_node_func(psrc);
			child->parent = cur;	
			stack_push(cur->children,child);
			stack_push((pStack)pdest->data,child);
		}
	}
  va_end(ap);
}

pTree tree_copy(pTree psrc, 
		        tree_func_node_new new_node_func, 
				tree_func_node_copy copy)
{
	pTree pdest = (pTree)memory_new(1, psrc->size);
	memcpy(pdest,psrc, sizeof(psrc->size));
	pdest->count = psrc->count;
	pdest->nlist = NULL;
	pdest->prefix = NULL;
	pdest->postfix = NULL;
	pdest->data = NULL;
	pdest->size = psrc->size;
	(*psrc->count)++;
	pdest->root = new_node_func(psrc);
	pdest->data = stack_new();
	stack_push( (pStack)pdest->data, pdest->root );
	tree_traveral_prefix(psrc, __tree_copy, pdest, new_node_func, copy);
	stack_free((pStack)pdest->data);
	pdest->data = psrc->data; 
	return pdest;
}

void __tree_build_node_list(pTree ptree, pTreeNode ptnode, va_list ap1)
{
	arraylist_add((pArrayList)ptree->nlist,ptnode);
}

void __tree_build_prefix_node_list(pTree ptree, pTreeNode ptnode, va_list ap1)
{
	arraylist_add((pArrayList)ptree->prefix,ptnode);
}

void __tree_build_postfix_node_list(pTree ptree, pTreeNode ptnode, va_list ap1)
{
	arraylist_add((pArrayList)ptree->postfix,ptnode);
}

void tree_build_node_list(pTree ptree)
{
	if ( ptree->nlist )
	{
		arraylist_free(ptree->nlist, NULL);			
		arraylist_free(ptree->postfix, NULL);			
		arraylist_free(ptree->prefix, NULL);			
	}
	// Even index ( 0, 2, 4, .... ) : leaf
	ptree->nlist = arraylist_new(100);
	ptree->postfix = arraylist_new(100);
	ptree->prefix = arraylist_new(100);

	tree_traveral_infix(ptree, __tree_build_node_list);
	tree_traveral_postfix(ptree, __tree_build_postfix_node_list);
	tree_traveral_prefix(ptree, __tree_build_prefix_node_list);

//	arraylist_trim((pArrayList)ptree->nlist);
	int i = 0;
	for ( i = 0 ; i < ptree->nlist->size ; i++ )
	{
		((pTreeNode)ptree->nlist->array[i])->id = i;
	}
}
