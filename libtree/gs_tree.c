#include<tree.h>


pTreeNode gene_tree_new_empty_node(pTree ptree)
{
	pGeneTreeNode pgnode = (pGeneTreeNode)memory_new(1, sizeof(GeneTreeNode));
	phylogeny_clear_node((pPhylogenyNode)pgnode);
	pgnode->vsubspecies = NULL;
	pgnode->psnode = NULL;
	return (pTreeNode)pgnode;
}

pTreeNode species_tree_new_empty_node(pTree ptree)
{
	pSpeciesTreeNode psnode = (pSpeciesTreeNode)memory_new(1, sizeof(SpeciesTreeNode));
	phylogeny_clear_node((pPhylogenyNode)psnode);
	psnode->vsubspecies = NULL;
	psnode->mask = NULL;
	psnode->gain = 0;
	psnode->numduplicated = 0;
	psnode->size = -1;
	return (pTreeNode)psnode;
}

pTree __gene_tree_new(tree_func_node_new new_tree_node_func, va_list ap1)
{
	pGeneTree pgtree = (pGeneTree) memory_new(1,sizeof(GeneTree));
	pgtree->super.size = sizeof(GeneTree);
	pgtree->species_tree = NULL;
	tree_new_fill((pTree)pgtree, gene_tree_new_empty_node );
	return (pTree)pgtree;
}

void gene_tree_node_copy(pTreeNode psrc, pTreeNode pdest )
{
	phylogeny_node_copy(psrc,pdest);
	((pGeneTreeNode)pdest)->vsubspecies = ((pGeneTreeNode)psrc)->vsubspecies;
	((pGeneTreeNode)pdest)->psnode = ((pGeneTreeNode)psrc)->psnode;
}

void species_tree_node_copy(pTreeNode psrc, pTreeNode pdest )
{
	phylogeny_node_copy(psrc,pdest);
	pSpeciesTreeNode psnode = (pSpeciesTreeNode)pdest;
	psnode->vsubspecies = ((pSpeciesTreeNode)psrc)->vsubspecies;
	psnode->gain = 0;
	psnode->mask = ((pSpeciesTreeNode)psrc)->mask;
}

extern void __phylogeny_free_node(pTree ptree, pTreeNode ptnode, va_list ap1);

void __gene_tree_free_node(pTree ptree, pTreeNode ptnode, va_list ap1)
{
  va_list ap;
  if (ap1) va_copy(ap, ap1);
	if ( *ptree->count == 0  )
	{
		pSpeciesTreeNode pnode = (pSpeciesTreeNode)ptnode;
		if ( pnode->vsubspecies ) vector_free(pnode->vsubspecies,NULL);
	}
	__phylogeny_free_node(ptree,ptnode,ap);		
  if (ap1) va_end(ap);
}

void __species_tree_free_node(pTree ptree, pTreeNode ptnode, va_list ap1)
{
  va_list ap;
  if (ap1) va_copy(ap, ap1);
	pSpeciesTreeNode pnode = (pSpeciesTreeNode)ptnode;
	if ( *ptree->count == 0  )
	{
		if ( pnode->vsubspecies ) vector_free(pnode->vsubspecies,NULL);
	}
	memory_free( pnode->mask );
	pnode->mask = NULL;
	__phylogeny_free_node(ptree,ptnode,ap);		
  if (ap1) va_end(ap);
}

void gene_tree_free(pGeneTree pgtree)
{
	int i;
	pArrayList nlist = pgtree->super.nlist;

	for ( i = 0 ; i < nlist->size ; i++ )
	{
		__gene_tree_free_node((pTree)pgtree, (pTreeNode)nlist->array[i], NULL );
	}
	memory_free(pgtree->species_tid);
	pgtree->species_tid = NULL;
	tree_free((pTree)pgtree);
}

void species_tree_free(pTree ptree)
{
	int i;
	pArrayList nlist = ptree->nlist;
	for ( i = 0 ; i < nlist->size ; i++ )
	{
		__species_tree_free_node((pTree)ptree, (pTreeNode)nlist->array[i], NULL );
	}
	tree_free(ptree);
}


size_t __taxaid_cmp(const void* arg1, const void* arg2)
{
	return (uintptr_t)arg1 - (uintptr_t)arg2;	
}

void __gene_species_tree_make_subspecies(pTree ptree)
{
	int i;

	if ( ptree->nlist->size == 2 ) return;	

	for ( i = 0 ; i < ptree->nlist->size ; i+=2 )
	{
		pSpeciesTreeNode psnode = (pSpeciesTreeNode)ptree->nlist->array[i];
		uintptr_t taxaid = psnode->super.taxaid;
		psnode = (pSpeciesTreeNode)((pTreeNode)psnode)->parent;
		while( psnode )
		{
			if ( psnode->vsubspecies == NULL )
			{
				psnode->vsubspecies = vector_new();
			}
			pVector subs = psnode->vsubspecies;
			uintptr_t tid = (uintptr_t)vector_get_by_cmp( subs, (void*)taxaid, __taxaid_cmp  );
			if ( tid )  break;
			vector_add(subs,(void*)taxaid);
			psnode = (pSpeciesTreeNode)((pTreeNode)psnode)->parent;
		}
	}
}

int gene_tree_map_species(pGeneTree pgtree)
{
	pArrayList gnlist = pgtree->super.nlist;
	pArrayList snlist = pgtree->species_tree->nlist;
	tree_clear_reg((pTree)pgtree);
	int i,j;

	if ( gnlist->size == 2 )
	{
		pGeneTreeNode pgnode = (pGeneTreeNode)gnlist->array[0];
		pSpeciesTreeNode psnode = NULL;
		int tid = (int)pgnode->super.taxaid;
		for ( j = 0 ; j < snlist->size ; j+=2 )
		{
			pPhylogenyNode p = (pPhylogenyNode)snlist->array[j];
			if ( p->taxaid  == tid ) 
			{
				psnode = (pSpeciesTreeNode)p;
				break;
			}
		}
		if ( psnode == NULL ) 
		{
			fprintf(stderr, "Error: There is no such the species( taxa id %d ) in the species tree\n",  tid );
			return -1;
		}
		pgnode->psnode = psnode;
		return 0;
	}

	// For leaf mapping
	for ( i = 0 ; i < gnlist->size ; i+=2 )
	{
		pGeneTreeNode pgnode = (pGeneTreeNode)gnlist->array[i];
		int tid = (int)pgnode->super.taxaid;
		for ( j = 0 ; j < snlist->size ; j+=2 )
		{
			pSpeciesTreeNode p = (pSpeciesTreeNode)snlist->array[j];
			if ( p->super.taxaid  == tid ) 
			{
				pgnode->psnode = p;
				break;
			}
		}
	}

	// For node mapping
	for ( i = 1 ; i < gnlist->size ; i+=2 )
	{
		pGeneTreeNode pgnode = (pGeneTreeNode)gnlist->array[i];
		pSpeciesTreeNode psnode = NULL;
		uintptr_t tid = (uintptr_t)pgnode->vsubspecies->head->data;
		// Find start point from species tree.
		for ( j = 0 ; j < snlist->size ; j+=2 )
		{
			pPhylogenyNode p = (pPhylogenyNode)snlist->array[j];
			if ( p->taxaid  == tid ) 
			{
				psnode = (pSpeciesTreeNode)p;
				break;
			}
		}
		if ( psnode == NULL ) 
		{
			fprintf(stderr, "Error: There is no such the species( taxa id %ld ) in the species tree\n",  tid );
			return -1;
		}
		if ( pgnode->vsubspecies->size > 1 )
		{
			int cnt = 1;
			psnode = (pSpeciesTreeNode)((pTreeNode)psnode)->parent;
			pLinkedListItem cur = pgnode->vsubspecies->head->next;
			while( cnt < pgnode->vsubspecies->size )
			{
				if ((uintptr_t)vector_get_by_cmp( psnode->vsubspecies, cur->data, __taxaid_cmp  ) )
				{
					cnt++;
					cur = cur->next;
				}
				else
				{
					psnode = (pSpeciesTreeNode)((pTreeNode)psnode)->parent;
				}
			}
		}
		pgnode->psnode = psnode;
	}
	return 0; 
}

int __gene_tree_get_index_of_taxon(pGeneTree pgtree, int tid, int numspecies )
{
	int k;
	for ( k = 0 ; k < numspecies ; k++ )
	{
		if ( pgtree->species_tid[k] == tid ) return k;
	}
	fprintf(stderr,"Error: There is no such taxa id, %d\n", tid);
	return -1;
}
	
pTree species_tree_new(char* sztree, phylogeny_func_parse_node fparse )
{
	pTree ptree = phylogeny_load_from_string(sztree, tree_new, species_tree_new_empty_node, fparse );
	if ( ptree == NULL ) return NULL;
	tree_build_node_list(ptree);
	__gene_species_tree_make_subspecies(ptree);

	int i;
	int numspecies = (ptree->nlist->size+1)/2;
	for ( i = 1 ;  i < ptree->nlist->size ; i+=2 )
	{
		pSpeciesTreeNode psnode = (pSpeciesTreeNode)ptree->nlist->array[i];		
		psnode->mask = (int*) memory_new( numspecies, sizeof(int) );
		//memset( psnode->mask, 0, sizeof(int)*numspecies);
	}
	for ( i = 0 ; i< ptree->nlist->size; i+=2 )
	{
		pSpeciesTreeNode psnode = (pSpeciesTreeNode)ptree->nlist->array[i];
		psnode = (pSpeciesTreeNode)((pTreeNode)psnode)->parent;
		while( psnode )
		{
			psnode->mask[i/2] = 1;
			psnode = (pSpeciesTreeNode)((pTreeNode)psnode)->parent;
		}
	}
	return ptree;
}

void __with_reg(pString pstr, pPhylogenyNode pnode)
{
	if ( pnode->name ) string_fadd(pstr,"%s", pnode->name );
	string_fadd(pstr,"<%d>", pnode->super.reg );
}

pGeneTree gene_tree_new(char* sztree, pTree pstree, phylogeny_func_parse_node fparse )
{
	pGeneTree ptree = (pGeneTree)phylogeny_load_from_string(sztree, __gene_tree_new, gene_tree_new_empty_node, fparse );
	if ( ptree == NULL ) return NULL;
	tree_build_node_list((pTree)ptree);
	ptree->species_tree = pstree;

	int i, j;
	if ( !pstree ) return ptree;
	int num_species = (pstree->nlist->size+1)/2;
	int* staxaid = (int*) memory_new(num_species, sizeof(int));
	for ( i = 0, j = 0 ; i < pstree->nlist->size ; i+=2 )
	{
		staxaid[j++] = ((pPhylogenyNode)pstree->nlist->array[i])->taxaid;
	}
	pArrayList nlist = ptree->super.nlist;
	for ( i = 0 ; i < nlist->size ; i+= 2)
	{
		pPhylogenyNode pgnode = (pPhylogenyNode)nlist->array[i];
		for( j = 0 ; j < num_species ; j++ )
		{
			if ( staxaid[j] == pgnode->taxaid )
			{
				pTreeNode p = (pTreeNode)pgnode;
				while( p && p->reg == 0)
				{
					p->reg = 1;
					p = p->parent;
				}
				break;
			}
		}
	}
	ptree->species_tid = staxaid;
	if ( phylogeny_delete_nodes_by_reg((pTree)ptree, __gene_tree_free_node ) == -1 )
	{
		return NULL;
	}
	__gene_species_tree_make_subspecies((pTree)ptree);
	return ptree;
}

void __species_tree_string_with_subspecies(pString pstr, pPhylogenyNode ppnode )
{
	pSpeciesTreeNode pnode = (pSpeciesTreeNode)ppnode;
	if ( ppnode->name ) string_add(pstr, ppnode->name );
	if ( pnode->vsubspecies && pnode->vsubspecies->size )
	{
		int i;
		pArrayList pal = vector_to_arraylist( pnode->vsubspecies );
		string_fadd(pstr,"<%ld", (uintptr_t)pal->array[0] );
		for ( i = 1 ; i < pal->size ; i++ )
		{
			string_fadd(pstr,",%ld", (uintptr_t)pal->array[i] );
		}
		string_add(pstr,">");
		arraylist_free(pal, NULL);
	}
}

void __gene_tree_string_with_subspecies(pString pstr, pPhylogenyNode ppnode )
{
	pSpeciesTreeNode pnode = ((pGeneTreeNode)ppnode)->psnode;
	if ( ppnode->name ) string_add(pstr, ppnode->name );
	if ( pnode && pnode->vsubspecies && pnode->vsubspecies->size )
	{
		int i;
		pArrayList pal = vector_to_arraylist( pnode->vsubspecies );
		string_fadd(pstr,"<%ld", (uintptr_t)pal->array[0] );
		for ( i = 1 ; i < pal->size ; i++ )
		{
			string_fadd(pstr,",%ld", (uintptr_t)pal->array[i] );
		}
		string_add(pstr,">");
		arraylist_free(pal, NULL);
	}
	else if ( pnode )
	{
		string_fadd(pstr,"<%d>", pnode->super.taxaid );
	}
}

pString gene_tree_string(pGeneTree ptree, int mode)
{
	return mode == TREE_DEBUG ?
		phylogeny_string_nhx((pTree)ptree, __gene_tree_string_with_subspecies, 0 ):
		phylogeny_string_nhx((pTree)ptree, NULL, 0 );
}

pString species_tree_string(pTree ptree, int mode)
{
	return mode == TREE_DEBUG ?
		phylogeny_string_nhx(ptree, __species_tree_string_with_subspecies, 0 ) :
		phylogeny_string_nhx(ptree, NULL, 0 );
}

void __gene_tree_build_sub_nodes_list(pTree ptree, pTreeNode ptnode, va_list ap1)
{
  va_list ap;
  va_copy(ap, ap1);
	pArrayList pal = va_arg(ap, pArrayList);	
	arraylist_add(pal,ptnode);
  va_end(ap);
}

void __gene_tree_fill_mark( pGeneTree pgtree, pVector vss, int numspecies, int* mark )
{
	uintptr_t tid;
	vector_rewind(vss);
	while( (tid=(uintptr_t)vector_next(vss)) )
	{
		tid = __gene_tree_get_index_of_taxon(pgtree, (int)tid, numspecies );
		mark[tid] = 1;		
	}
}

void __gene_tree_duplication_count_block_sub_count(pGeneTree pgtree, pGeneTreeNode pgnode, pArrayList nsublist, int* mark, int numspecies )
{
	int k;
	Tree tmp; tmp.nlist=  NULL;
	tmp.root = (pTreeNode)pgnode;
	arraylist_clear( nsublist );
	tree_traveral_infix(&tmp, __gene_tree_build_sub_nodes_list, nsublist );

	for ( k = 1 ; k < nsublist->size ; k+=2 )
	{
		pGeneTreeNode pg = (pGeneTreeNode)nsublist->array[k];
		if ( pg->super.duplicated == PHYLOGENY_SPECIATION_DUPLICATED && pg->vsubspecies->size > 1 )
		{
			__gene_tree_fill_mark(pgtree, pg->psnode->vsubspecies, numspecies, mark );
		}
	}
}

void __gene_tree_duplication_count_gain(pTree pstree, pSpeciesTreeNode proot, pArrayList nsublist, int* mark, int numspecies)
{
	int k;
	Tree tmp; tmp.nlist=  NULL;
	tmp.root = (pTreeNode)proot;
	arraylist_clear( nsublist );
	tree_traveral_postfix(&tmp, __gene_tree_build_sub_nodes_list, nsublist );

	for ( k = 0 ; k < numspecies ; k++ )
	{
		if( proot->mask[k] && mark[k] == 0  )
		{
			pTreeNode ps = (pTreeNode)pstree->nlist->array[k*2];
			ps->reg = 1;
			ps->parent->reg++;
		}
	}

	for ( k = 0 ; k < nsublist->size ; k++ )
	{
		pTreeNode pg = (pTreeNode)nsublist->array[k];
		if ( tree_is_leaf(pg) ) continue;
		if ( pg->reg == 2 ) 
		{
			pg->reg = 1;
			((pTreeNode)pg->children->head->data)->reg = 0;
			((pTreeNode)pg->children->tail->data)->reg = 0;
			if ( pg->parent ) pg->parent->reg++;
		}
		else
		{
			pg->reg = 0;
		}
	}
	for ( k = 0 ; k < nsublist->size ; k++ )
	{
		pSpeciesTreeNode ps = (pSpeciesTreeNode)nsublist->array[k];
		if ( ((pTreeNode)ps)->reg ) ps->gain--;
	}
	
	for ( k = 0 ; k < nsublist->size ; k++ )
	{
		pTreeNode ps = (pTreeNode)nsublist->array[k];
		ps->reg = 0;
	}
}

void species_tree_clear_gain(pTree pstree)
{
	int i;
	pArrayList nslist = pstree->nlist;
	for ( i = 0 ; i < nslist->size ; i++ )
	{
		pSpeciesTreeNode psnode = (pSpeciesTreeNode)nslist->array[i];
		{
			psnode->gain = 0;
			psnode->numduplicated = 0;
		}
	}
}

int gene_tree_duplication_count(pGeneTree pgtree)
{
	int i,k;
	pTree pstree = pgtree->species_tree;
	pArrayList nslist = pstree->nlist;
	pArrayList nglist = pgtree->super.nlist;

	int numspecies = (nslist->size+1)/2;

	int* mark = (int*) memory_new( numspecies, sizeof(int) );
	pArrayList nsublist = arraylist_new(100);

	if ( nglist->size == 2 )
	{
		pGeneTreeNode pgnode = (pGeneTreeNode)nglist->array[0];
		memset( mark, 0, sizeof(int) * numspecies );
		int k = __gene_tree_get_index_of_taxon(pgtree, pgnode->super.taxaid, numspecies );
		mark[k] = 1;		
		__gene_tree_duplication_count_gain(pstree,(pSpeciesTreeNode)pstree->root,nsublist,mark,numspecies);
		memory_free(mark);
		mark = NULL;
		arraylist_free( nsublist, NULL );
		return 0;
	}

	for ( i = 0 ; i < nglist->size ; i++ )
	{
		pGeneTreeNode pgnode = (pGeneTreeNode)nglist->array[i];
		pSpeciesTreeNode psnode = pgnode->psnode;
		pGeneTreeNode pgparent = (pGeneTreeNode)((pTreeNode)pgnode)->parent;
		pSpeciesTreeNode proot = pgparent ? pgparent->psnode : psnode;

		if ( pgnode->super.duplicated == PHYLOGENY_SPECIATION_DUPLICATED ) 
		{
			psnode->gain++;
			psnode->numduplicated++;
			if ( ( pgparent && pgparent->super.duplicated == PHYLOGENY_SPECIATION_DUPLICATED ) )
			{
				if ( pgparent->psnode == psnode ) continue;
				memset( mark, 0, sizeof(int) * numspecies );
				__gene_tree_fill_mark( pgtree, pgnode->vsubspecies, numspecies, mark );
				__gene_tree_duplication_count_block_sub_count(pgtree,pgnode,nsublist,mark,numspecies);
				__gene_tree_duplication_count_gain(pstree,proot,nsublist,mark,numspecies);
			}
			continue;
		}
		if ( pgparent && pgparent->super.duplicated != PHYLOGENY_SPECIATION_DUPLICATED ) continue;
		memset( mark, 0, sizeof(int) * numspecies );
		if ( tree_is_leaf((pTreeNode)pgnode))
		{
			if ( pgparent->vsubspecies->size == 1 ) continue;
			k = __gene_tree_get_index_of_taxon(pgtree, pgnode->super.taxaid, numspecies );
			mark[k] = 1;		
		}
		else
		{
			__gene_tree_fill_mark( pgtree, pgnode->vsubspecies, numspecies, mark );
			__gene_tree_duplication_count_block_sub_count(pgtree,pgnode,nsublist,mark,numspecies);
		}
		__gene_tree_duplication_count_gain(pstree,proot,nsublist,mark,numspecies);
	}

	pGeneTreeNode pgroot = (pGeneTreeNode)((pTree)pgtree)->root;
	if ( pgroot->super.duplicated != PHYLOGENY_SPECIATION_DUPLICATED &&
		 pgroot->psnode != (pSpeciesTreeNode)pstree->root ) 
	{
		memset( mark, 0, sizeof(int) * numspecies );
		__gene_tree_fill_mark( pgtree, pgroot->vsubspecies, numspecies, mark );
		__gene_tree_duplication_count_block_sub_count(pgtree,pgroot,nsublist,mark,numspecies);
		__gene_tree_duplication_count_gain(pstree,(pSpeciesTreeNode)pstree->root,nsublist,mark,numspecies);
	}

	memory_free(mark);
	mark = NULL;
	arraylist_free(nsublist, NULL);
	return 0;
}

void __species_tree_with_id(pString pstr, pPhylogenyNode pnode)
{
	pSpeciesTreeNode psnode = (pSpeciesTreeNode)pnode;
	if ( pnode->name )
	{
		string_fadd(pstr,"%s", pnode->name);
		if ( psnode->size != -1 ) string_fadd(pstr,"_%d", psnode->size);
		string_fadd(pstr,"<%d/%d>", psnode->numduplicated,
				psnode->numduplicated - psnode->gain );
	}
	else
	{
		string_fadd(pstr,"<%d/%d>", 
				psnode->numduplicated,
				psnode->numduplicated - psnode->gain );
	}
}

void gene_tree_duplication_count_report(pTree pstree, FILE* fout)
{
	pString psp = phylogeny_string_newick(pstree, __species_tree_with_id,PS_SKIP_BL);
	fprintf( fout, "%s\n", psp->buf );
	string_free(psp);
}

double gstree_mp_remark(pString pstr, pTree ptree, pMetapostConfig pmc, va_list ap1 )
{
  va_list ap;
  va_copy(ap, ap1);
	char* title = va_arg(ap, char*);
	string_fadd( pstr, "label( btex  etex, (0.1u, %fu));\n", pmc->height + 0.7 ) ;
	string_fadd( pstr, "label( btex %s etex, (0.1u, %fu));\n", title, pmc->height + 0.5 ) ;
  va_end(ap);
	return 0;	
}

double gene_tree_mp_annotation(pString pstr, pTreeNode pnode, pMetapostConfig pmc, va_list ap1 )
{
    string_add( pstr, ";\n");
	if ( ((pPhylogenyNode)pnode)->branchlength == -1 ) return 0;
	if ( pnode->parent )
	{
		string_fadd( pstr, "ypart mid[%d] = ypart(p[%d]);\n", pnode->id, pnode->id );
		string_fadd( pstr, "xpart mid[%d] = (xpart(p[%d])+xpart(p[%d]))/2;\n", pnode->id, pnode->id, pnode->parent->id );
		string_fadd( pstr, "label.top( btex $%g$ ", ((pPhylogenyNode)pnode)->branchlength );
		string_fadd( pstr, "etex, mid[%d]);\n", pnode->id  );
	}
	return 0;
}

double species_tree_mp_annotation(pString pstr, pTreeNode pnode, pMetapostConfig pmc, va_list ap1 )
{
	pSpeciesTreeNode psnode = (pSpeciesTreeNode)pnode;
	string_add(pstr, ";\n");
	double last = 0;
	if ( psnode->size >= 0 )
	{
		string_fadd( pstr, "label.urt( btex	%d etex, p[%d]);\n", psnode->size, pnode->id );
	}
	if ( pnode->parent )
	{
		string_fadd( pstr, "xpart mid[%d] = xpart(p[%d]);\n", pnode->id, pnode->id );
        string_fadd( pstr, "ypart mid[%d] = (ypart(p[%d])+ypart(p[%d]))/2;\n", pnode->id, pnode->id, pnode->parent->id );
		string_fadd( pstr, "label.rt( btex gain = %d ", psnode->numduplicated  );
		string_fadd( pstr, "etex, mid[%d]);\n", pnode->id  );
		string_fadd( pstr, "label.rt( btex loss = %d ", psnode->numduplicated - psnode->gain );
		last -= 0.15;
		string_fadd( pstr, "etex, mid[%d] + (0,%fu));\n",  pnode->id, last );
	}
	else
	{
		string_fadd( pstr, "label.rt( btex gain = %d ", psnode->numduplicated  );
		string_fadd( pstr, "etex, p[%d] + (0,%fu));\n", pnode->id, 0.15 );
	}
	return last;
}

pString gene_tree_metapost(pGeneTree pgtree, int id, char* title, double width, double height )
{
	MetapostConfig mc;
	mc.id = id;
	mc.unit = MP_UNIT_IN;
	mc.dir = MP_DIR_HORIZONTAL;
	mc.shape = MP_SHAPE_RECT | MP_SHAPE_MOST_CENTER;
	mc.fmod = gene_tree_mp_annotation;
	mc.fremark = gstree_mp_remark;
	mc.width = width;
	mc.height = height;
	return phylogeny_to_mp( (pTree)pgtree, &mc, title );
}


pString species_tree_metapost(pTree pstree, int id, char* title, double width, double height )
{
	MetapostConfig mc;
	mc.id = id;
	mc.unit = MP_UNIT_IN;
	mc.dir = MP_DIR_VERTICAL;
	mc.shape = MP_SHAPE_RECT | MP_SHAPE_MOST_CENTER;
	mc.fmod = species_tree_mp_annotation;
	mc.fremark = gstree_mp_remark;
	mc.width = width;
	mc.height = height;
	return phylogeny_to_mp( (pTree)pstree, &mc, title );
}



