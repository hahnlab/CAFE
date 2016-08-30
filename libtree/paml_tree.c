#include<tree.h>
#include<io.h>

extern void __phylogeny_free_node(pTree ptree, pTreeNode ptnode, va_list ap1);
pTreeNode paml_tree_new_empty_node(pTree ptree)
{
	pPamlNode pnode = (pPamlNode)memory_new(1,sizeof(PamlNode));	
	pnode->N = -1;
  pnode->accudS = -1;
	phylogeny_clear_node((pPhylogenyNode)pnode);
	return (pTreeNode)pnode;
}

pTree __paml_tree_new(tree_func_node_new new_tree_node_func, va_list ap1)
{
	pPamlTree ptree = (pPamlTree) memory_new(1, sizeof(PamlTree));
	tree_new_fill((pTree)ptree, paml_tree_new_empty_node);
	ptree->super.size = sizeof(PamlTree);
	return (pTree)ptree;
}

pPamlTree paml_tree_new(char* file)
{
	pPamlTree ptree = NULL;
	FILE* fp;
	if (  ( fp = fopen(file,"r") ) == NULL )
	{
		fprintf(stderr, "ERROR(paml): Cannot open %s\n", file );
		perror("Error\n");
		return NULL;
	}

	pString pstr = string_new();

	int num_genes = -1;
	char** gene_names = NULL;
	double lnL = 0;
	int genecnt = 0;
	double** nei_gojobori_w = NULL;
  double** nei_gojobori_dN = NULL;
  double** nei_gojobori_dS = NULL;
  
	while( file_read_line(pstr, fp) > 0 )
	{
		if ( num_genes == -1 && strncmp( pstr->buf, "ns", 2 ) == 0 )
		{
			sscanf( pstr->buf, "ns = %d", &num_genes );	
			gene_names = (char**) memory_new(num_genes, sizeof(char*));
		}
		else if ( pstr->buf[0] == '#' )
		{
      if ( strncmp( pstr->buf, "# site patterns", 15 ) == 0 ) {
        // version 3.13 output, ignore
        continue;
      }
			char* n = index(pstr->buf, ' ');
			*index(++n,' ') = '\0';
			gene_names[genecnt] = (char*) memory_new(strlen(n)+1, sizeof(char));
			strcpy( gene_names[genecnt], n );
			genecnt++;
		}
		else if ( strncmp( pstr->buf, "lnL", 3) == 0 )
		{
			char* n = index( pstr->buf, '-' );
			*index(n,' ') = '\0';
			sscanf( n, "%lf", &lnL );
		}
		else if ( strncmp(pstr->buf, "Nei &", 5 ) == 0 )
		{
			int i,j;
			while(1)
			{
				file_read_line(pstr, fp);
				if ( pstr->buf[0] == '\n' ) break;
			}
			file_read_line(pstr, fp);
			//int num_ws = (int)(genecnt * ( genecnt - 1 ) / 2.0);
      
      nei_gojobori_w = calloc( genecnt, sizeof(*nei_gojobori_w) );
      nei_gojobori_w[0] = calloc( genecnt*genecnt, sizeof(**nei_gojobori_w) );
      for ( i=1; i<genecnt; i++) {
        nei_gojobori_w[i] = nei_gojobori_w[0] + i * genecnt;
      }
      nei_gojobori_dN = calloc( genecnt, sizeof(*nei_gojobori_dN) );
      nei_gojobori_dN[0] = calloc( genecnt*genecnt, sizeof(**nei_gojobori_dN) );
      for ( i=1; i<genecnt; i++) {
        nei_gojobori_dN[i] = nei_gojobori_dN[0] + i * genecnt;
      }
      nei_gojobori_dS = calloc( genecnt, sizeof(*nei_gojobori_dS) );
      nei_gojobori_dS[0] = calloc( genecnt*genecnt, sizeof(**nei_gojobori_dS) );
      for ( i=1; i<genecnt; i++) {
        nei_gojobori_dS[i] = nei_gojobori_dS[0] + i * genecnt;
      }  
			int sum = 0;
			for ( i = 1 ; i < genecnt; i++ )
			{
				file_read_line(pstr, fp);
				char* ptr = &pstr->buf[20];
				double w, dN, dS;
				for ( j = 0 ; j < i ; j++, sum++ )
				{
				 	sscanf( ptr, "%lf (%lf %lf)", &w, &dN, &dS );
					nei_gojobori_w[i][j] = w;
					nei_gojobori_dN[i][j] = dN;
					nei_gojobori_dS[i][j] = dS;
					ptr = index( ptr, ')');
					if ( ptr ) ptr++;
				}	
			}
		}
		else if ( ptree == NULL && strncmp( pstr->buf, "tree len", 8 ) == 0 )
		{
			file_read_line(pstr, fp);
			file_read_line(pstr, fp);
			ptree = (pPamlTree)phylogeny_load_from_string(pstr->buf,__paml_tree_new, paml_tree_new_empty_node, NULL );
			if ( ptree == NULL ) return NULL;
			// check if tree is unrooted and root the tree
			if (((pTree)ptree)->root->children->size > 2) {
				ptree = (pPamlTree)phylogeny_root_tree((pTree)ptree, paml_tree_new_empty_node, NULL );
			}
			tree_build_node_list((pTree)ptree);
			ptree->num_genes = num_genes;
			ptree->gene_names = gene_names;
			ptree->lnL = lnL;
			ptree->nei_gojobori_w = nei_gojobori_w;
			ptree->nei_gojobori_dN = nei_gojobori_dN;
			ptree->nei_gojobori_dS = nei_gojobori_dS;
			pArrayList nlist = ptree->super.prefix;
			ptree->nodelist = (pPamlNode*) memory_new( nlist->size, sizeof(pPamlNode) );
			int id = num_genes;
			int i;
			for ( i = 0 ; i < nlist->size; i++ )
			{
				pPhylogenyNode pnode = (pPhylogenyNode)nlist->array[i];
				if ( pnode->super.children )
				{
					ptree->nodelist[id++] = (pPamlNode)pnode;
				}
				else
				{
					int idx = atoi(pnode->name)-1;
					ptree->nodelist[idx] = (pPamlNode)pnode;
					memory_free(pnode->name);
					pnode->name = NULL;
					pnode->name = gene_names[idx];
				}
			}
		}
		else if ( strncmp(pstr->buf, "dN/dS", 5) == 0 )
		{
			char* n = index(pstr->buf, '=' );
			*index(++n,')') = '\0';
			ptree->model = atoi(n);
			if( ptree->model == 4 )
			{
				int i;
				ptree->omega = (double*) memory_new( ptree->model*2, sizeof(double) );
				ptree->prob = (double*) memory_new( ptree->model, sizeof(double) );
				while(1)
				{
					file_read_line(pstr, fp);
					if ( strncmp(pstr->buf, "prop", 4 ) == 0 )
					{
						pArrayList list = string_pchar_space_split(pstr->buf);
						for( i = 1 ; i <= 4 ; i++ )
						{
							sscanf( (char*)list->array[i], "%lf", &ptree->prob[i-1] );
						}
						arraylist_free(list, NULL);
						file_read_line(pstr, fp);
						list = string_pchar_space_split(pstr->buf);
						for( i = 2 ; i <= 5 ; i++ )
						{
							sscanf( (char*)list->array[i], "%lf", &ptree->omega[i-2] );
						}
						arraylist_free(list, NULL);
						file_read_line(pstr, fp);
						list = string_pchar_space_split(pstr->buf);
						for( i = 2 ; i <= 5 ; i++ )
						{
							sscanf( (char*)list->array[i], "%lf", &ptree->omega[4+i-2] );
						}
						arraylist_free(list, NULL);
						break;
					}
				}
			}
			else
			{
				ptree->omega = (double*) memory_new( ptree->model, sizeof(double) );
				ptree->prob = (double*) memory_new( ptree->model, sizeof(double) );
				file_read_line(pstr, fp);
	
				file_read_line(pstr, fp);
				pArrayList list = string_pchar_space_split(pstr->buf);
				int i;
				for ( i = 1 ; i < list->size ; i++ )
				{
					ptree->prob[i-1] = atof((char*)list->array[i]);
				}
				arraylist_free(list, NULL);
	
				file_read_line(pstr, fp);

				int j=0;
        for ( i = 3 , j = 0; i < pstr->length && j < ptree->model; i+=9 )
        {       
          char t = pstr->buf[i+9];
          pstr->buf[i+9] = '\0'; 
          ptree->omega[j++] = atof(&pstr->buf[i]);
          pstr->buf[i+9] = t;
        }       
        
			}
		}
		else if ( strncmp(pstr->buf, " branch", 6 ) == 0 )
		{
			file_read_line(pstr, fp);
			do
			{
				file_read_line(pstr, fp);
				if( strlen(pstr->buf) < 2 ) break;
				pstr->buf[64] = '\0';
				pArrayList list = string_pchar_space_split(pstr->buf);
				pPamlNode pnode;
				int from, to;
				sscanf( (char*)list->array[0], "%d..%d", &from,&to);
				from--; to--;
				if ( from < num_genes ) pnode = ptree->nodelist[from];
				else if ( to < num_genes || from < to ) pnode = ptree->nodelist[to];
				else pnode = ptree->nodelist[from];
				pnode->N = atof((char*)list->array[2]);
				pnode->S = atof((char*)list->array[3]);
				pnode->omega = atof((char*)list->array[4]);
				pnode->dN = atof((char*)list->array[5]);
				pnode->dS = atof((char*)list->array[6]);
				ptree->meanOmega = pnode->omega;
			}while( 1 );
		}
	}
	string_free(pstr);
	fclose(fp);
	return ptree;
}


double paml_tree_mp_remark(pString pstr, pTree ptree, pMetapostConfig pmc, va_list ap1)
{
  va_list ap;
  va_copy(ap, ap1);
	pPamlTree pptree = (pPamlTree)ptree;
	char* title = va_arg(ap,char*);
	string_add( pstr, "\n% annotation\n");
	string_fadd( pstr, "label( btex %s : %s etex, (0.1u, %fu));\n", title,  (pptree->model == 2 ? "Netural" : "Positive"), pmc->height + 1 ) ;
	string_fadd( pstr, "label( btex $\\ln L$ : %f etex, (0.1u, %fu));\n", pptree->lnL, pmc->height + 0.8 ) ;
	string_fadd( pstr, "label( btex $\\bar{\\omega}$ : %f etex, (0.1u, %fu));\n", pptree->meanOmega, pmc->height + 0.60  ) ;
	static char* omega[] = { "\\omega < 1" , "\\omega = 1", "\\omega > 1" };

	int i;
	for ( i = 0 ; i < pptree->model ; i++ )
	{
		string_fadd( pstr, "label( btex $%s$ : %f(%3.2f) etex, (0.1u, %fu));\n", omega[i],pptree->omega[i], pptree->prob[i], pmc->height + 0.45 - i * 0.15  ) ;
	}
	

/*
	string_fadd( pstr, "label( btex \\begin{tabular}{|c|c|c|} \\hline \
				            & \\omega & P(\\omega) \\\\ \\hline \
				\\omega < 1 &  %f & %f%  \\\\ \
				\\omega == 1 &  %f & %f% \\\\ ", 
				pptree->omega[0], pptree->prob[0],
				pptree->omega[1], pptree->prob[1] );
	if ( pptree->model == 3 )
	{
		string_fadd( pstr, "\\omega > 1 &  %f & %f% \\\\ ", pptree->omega[2], pptree->prob[2] );
	}
	string_fadd( pstr ,"\\hline\n\\end{taular} etex, (0.1u, %fu));\n", title, pmc->height ) ;
*/
  va_end(ap);
	return 0;
}

double paml_tree_mp_annotation(pString pstr, pTreeNode ptnode, pMetapostConfig pmc, va_list ap1)
{
	pPamlNode pnode = (pPamlNode)ptnode; 

	double last = 0;
    string_add( pstr, ";\n");
	if ( ptnode->parent )
	{
		string_fadd( pstr, "ypart mid[%d] = ypart(p[%d]);\n", ptnode->id, ptnode->id );
		string_fadd( pstr, "xpart mid[%d] = (xpart(p[%d])+xpart(p[%d]))/2;\n", ptnode->id, ptnode->id, ptnode->parent->id );
		if ( pnode->super.branchlength > 0 )
		{
			string_fadd( pstr, "label.rt( btex $\\small{l = %g}$ etex, mid[%d] + (-0.2u,%fu));\n", pnode->super.branchlength, ptnode->id, 0.50 );
		}
		if ( pnode->N != -1 )
		{
			string_fadd( pstr, "label.rt( btex $\\small{dN = %g}$ etex, mid[%d] + (-0.2u,%fu));\n",  pnode->dN, ptnode->id, 0.35 );
			string_fadd( pstr, "label.rt( btex $\\small{dS = %g}$ etex, mid[%d] + (-0.2u,%fu));\n",  pnode->dS, ptnode->id, 0.20 );
			string_fadd( pstr, "label.rt( btex $\\small{accuS = %g}$ etex, mid[%d] + (-0.2u,%fu));\n",  pnode->accudS, ptnode->id, 0.05 );
		}
	}
	return last;
}

void paml_tree_free(pPamlTree ptree)
{
	if (ptree->nei_gojobori_w && ptree->nei_gojobori_w[0]) { memory_free(ptree->nei_gojobori_w[0]); ptree->nei_gojobori_w[0] = NULL; }
	if (ptree->nei_gojobori_dN && ptree->nei_gojobori_dN[0]) { memory_free(ptree->nei_gojobori_dN[0]); ptree->nei_gojobori_dN[0] = NULL;}
	if (ptree->nei_gojobori_dS && ptree->nei_gojobori_dS[0]) { memory_free(ptree->nei_gojobori_dS[0]); ptree->nei_gojobori_dS[0] = NULL;}
	if (ptree->nei_gojobori_w) { memory_free(ptree->nei_gojobori_w); ptree->nei_gojobori_w = NULL;}
	if (ptree->nei_gojobori_dN) { memory_free(ptree->nei_gojobori_dN); ptree->nei_gojobori_dN = NULL; } 
	if (ptree->nei_gojobori_dS) { memory_free(ptree->nei_gojobori_dS); ptree->nei_gojobori_dS = NULL; }
  if (ptree->gene_names) { memory_free(ptree->gene_names); ptree->gene_names = NULL;}
  if (ptree->omega) { memory_free(ptree->omega); ptree->omega = NULL;}
  if (ptree->prob) { memory_free(ptree->prob); ptree->prob = NULL;}
	if (ptree->nodelist) { memory_free(ptree->nodelist); ptree->nodelist = NULL;}
	phylogeny_free((pTree)ptree);
}

pString paml_tree_metapost(pPamlTree ptree, int id, char* title, double width, double height)
{
	pMetapostConfig pmc = (pMetapostConfig) memory_new(1,sizeof(MetapostConfig));
	pmc->id = id;
	pmc->unit = MP_UNIT_IN;
	pmc->dir = MP_DIR_HORIZONTAL;
	pmc->shape = MP_SHAPE_RECT | MP_SHAPE_MOST_CENTER;
	pmc->fmod = paml_tree_mp_annotation;
	pmc->fremark = paml_tree_mp_remark;
	pmc->width = width;
	pmc->height = height;
	pString pstr = phylogeny_to_mp( (pTree)ptree, pmc, title);
	memory_free(pmc);
	pmc = NULL;
	return pstr;
}

void __paml_tree_string(pString pstr, pPhylogenyNode ptnode)
{
	pPamlNode pnode = (pPamlNode)ptnode;
	if ( ptnode->name )
	{
		string_fadd( pstr, "%s", ptnode->name );
	}
	string_fadd( pstr, "<%f/%f>", pnode->dN, pnode->dS );
}

void paml_tree_print_info(pPamlTree ptree)
{
	switch( ptree->model )
	{
		case 2:
			fprintf(stderr, "Neutral:\n");
			break;
		case 3:
			fprintf(stderr, "Positive:\n");
			break;
	}		
	int i;
	for ( i = 0 ; i < ptree->model ; i++ )
	{
		fprintf(stderr, "%4.3f(%4.3f%%) ", ptree->omega[i], ptree->prob[i]*100 );
	}
	fprintf(stderr, "\n");
	pString pstr = paml_tree_string(ptree);
	fprintf(stderr, "%s", pstr->buf);
	string_free(pstr);
}

pString paml_tree_string(pPamlTree ptree)
{
	return phylogeny_string((pTree)ptree,__paml_tree_string);
}


void paml_tree_accu_dS_from_root(pPamlTree ptree)
{
  paml_tree_accu_dS_from_node(ptree, (pPamlNode)((pTree)ptree)->root, NULL);
}

void paml_tree_accu_dS_from_node(pPamlTree ptree, pPamlNode pstartnode, double** dS_matrix)
{
	int i = 0;
  int row = 0;
  //pTreeNode sibling = phylogeny_get_sibling((pTree)&(ptree->super), (pTreeNode)&(pstartnode->super.super));
	pPamlNode pnode = NULL;
	pArrayList prefix = ptree->super.prefix;
  // find the node to start accumulating dS
  // initialize starting node to 0
  do
  {
		pnode = (pPamlNode)prefix->array[i];          
    //if (pnode == (pPamlNode)sibling) {
      //break;
    //}
    pnode->accudS = 0;
    row = ((pTreeNode)pnode)->id;
    i++;
  }  while (pnode != pstartnode);
  // node found. start at i, accumulate dS
	for ( ; i < prefix->size; i++ )
	{
		pnode = (pPamlNode)prefix->array[i];   
    //if (pnode == (pPamlNode)sibling) {
    //  break;
    //}
		pPamlNode parent = (pPamlNode)(((pTreeNode)pnode)->parent);

    if (parent) 
    {
      if ( parent->accudS >= 0 ) 
      {
        pnode->accudS = parent->accudS+pnode->dS;
        if (dS_matrix) 
        {
          dS_matrix[row][((pTreeNode)pnode)->id] = pnode->accudS;
          dS_matrix[((pTreeNode)pnode)->id][row] = pnode->accudS;
        }
      }
      else 
      {
        fprintf(stderr, "ERROR(paml): parent accudS must be no less than zero\n" );
        perror("Error\n");
        return;      
      }      
    }
    else 
    {
      fprintf(stderr, "ERROR(paml): must have a parent\n" );
      perror("Error\n");
      return;      
    }
  }
  // now fill the blanks of the row using parent distance information.
  if (dS_matrix)
  {
    for (i=0; i< prefix->size; i++)
    {
      pnode = (pPamlNode)prefix->array[i];          
      if (dS_matrix[row][((pTreeNode)pnode)->id] < 0)
      {
        pPamlNode parent = (pPamlNode)(((pTreeNode)pnode)->parent);
        double a = dS_matrix[((pTreeNode)parent)->id][row];
        double b = dS_matrix[((pTreeNode)parent)->id][((pTreeNode)pnode)->id];
        if (a >=0 && b >=0) {
          dS_matrix[row][((pTreeNode)pnode)->id] = a + b;
          dS_matrix[((pTreeNode)pnode)->id][row] = a + b;
        }
        else {
          fprintf(stderr, "ERROR(paml): parent row should be full already\n");
          perror("Error\n");
          return;      
        }
      }
    }
  }
}

void paml_tree_accu_dS_from_bottom_with_longer_dS(pPamlTree ptree)
{
	pArrayList postfix = ptree->super.postfix;
	int i;
	for ( i = 0 ; i < postfix->size-1 ; i++ )
	{
		pPamlNode pnode = (pPamlNode)postfix->array[i];
		pPamlNode parent = (pPamlNode)(((pTreeNode)pnode)->parent);
    /* take larger value dS of two children */
    double accu;
    if (pnode->accudS < 0 ) 
    {
      pnode->accudS = 0;
    }
    accu = pnode->accudS + pnode->dS;
    if ( parent->accudS < accu )
		{
			parent->accudS = accu;
		}
  }
}

void paml_tree_accu_dS_from_bottom_with_average_dS(pPamlTree ptree)
{
	pArrayList postfix = ptree->super.postfix;
	int i;
	for ( i = 0 ; i < postfix->size-1 ; i++ )
	{
		pPamlNode pnode = (pPamlNode)postfix->array[i];
		pPamlNode parent = (pPamlNode)(((pTreeNode)pnode)->parent);
    /* take average value dS of two children */
    
    double accu;
    if (pnode->accudS < 0 ) 
    {
      accu = pnode->dS;
    }
    else 
    {
      accu = pnode->accudS + pnode->dS;
    }
    double parent_accu = parent->accudS;
    if (parent->accudS < 0) 
    {
      parent->accudS = accu;
    }
    else 
    {
      parent->accudS = ( parent_accu + accu )/2;
    }
  }
}

int paml_tree_get_nodes_over_threshold(pPamlTree ptree, double threshold, int* rtn)
{
	paml_tree_accu_dS_from_bottom_with_longer_dS(ptree);
	pArrayList nlist = ptree->super.nlist;
	int i;
	for( i = 1 ; i < nlist->size ; i+=2 )
	{
		pPamlNode pnode = (pPamlNode)nlist->array[i];
		if ( pnode->accudS < threshold )
		{
			((pTreeNode)pnode)->reg = 1;
		}
		else
		{
			((pTreeNode)pnode)->reg = 0;
		}
	}
	int cnt = 0;
	for( i = 1 ; i < nlist->size ; i+=2 )
	{
		pTreeNode pnode = (pTreeNode)nlist->array[i];
		if ( pnode->reg == 1 )
		{
			pTreeNode parent = pnode->parent;	
			while( parent && parent->reg == 1 )
			{
				pnode = parent;
				parent = pnode->parent;	
			}
			rtn[cnt++] = pnode->id;	
		}
	}
	return cnt;
}

int paml_tree_check_significant(pPamlTree ptree, pPamlTree parg, double t)
{
   return 2 * ( ptree->lnL - parg->lnL ) > t;
}

double paml_tree_average_w(pPamlTree ptree)
{
	int i;
	double avg = 0;
	for ( i = 0 ; i < ptree->model ; i++ )
	{
		avg += ptree->omega[i] * ptree->prob[i];
	}	
	return avg;
}

double paml_tree_average_nei_gojobori(pPamlTree ptree)
{
	int i, j;
	double sum = 0;
	if ( ptree->nei_gojobori_w == NULL ) return -1;
	//int cnt = (int)(ptree->num_genes * ( ptree->num_genes - 1 ) / 2.0);
	int realcnt = 0;
	for ( i = 1 ; i < ptree->num_genes ; i++ )
	{
    for ( j = 0 ; j < i ; j++ )
    {
      if ( ptree->nei_gojobori_w[i][j] != -1 )
      {
        sum += ptree->nei_gojobori_w[i][j];
        realcnt++;
      }
    }
	}
	return realcnt == 0 ? -1 : sum/realcnt; 
}

void __paml_tree_free_node(pTree ptree, pTreeNode ptnode, va_list ap)
{
  __phylogeny_free_node( ptree, ptnode, NULL);
}


void __paml_tree_delete_nodes(pPamlTree ptree, pPamlNode pnode, tree_func_node freenode ) 
{
  pTreeNode sibling = phylogeny_get_sibling((pTree)ptree,(pTreeNode)pnode);
	if ( tree_is_root((pTree)ptree, ((pTreeNode)pnode)->parent) ) 
	{
		freenode((pTree)ptree, ((pTree)ptree)->root, NULL);	
		((pTree)ptree)->root = sibling;	
		((pTree)ptree)->root->parent = NULL;	
    ((pPamlNode)((pTree)ptree)->root)->N = 0;
    ((pPamlNode)((pTree)ptree)->root)->S = 0;
    ((pPamlNode)((pTree)ptree)->root)->dS = 0;
    ((pPamlNode)((pTree)ptree)->root)->dN = 0;
    ((pPamlNode)((pTree)ptree)->root)->omega = 0;
    ((pPamlNode)((pTree)ptree)->root)->accudS = 0;
	}
	else
	{
		((pPhylogenyNode)sibling)->branchlength += ((pPhylogenyNode)sibling->parent)->branchlength;
    ((pPamlNode)sibling)->N += ((pPamlNode)sibling->parent)->N;
    ((pPamlNode)sibling)->S += ((pPamlNode)sibling->parent)->S;
    ((pPamlNode)sibling)->dN += ((pPamlNode)sibling->parent)->dN;
    ((pPamlNode)sibling)->dS += ((pPamlNode)sibling->parent)->dS;
    ((pPamlNode)sibling)->omega = ((pPamlNode)sibling)->dN/((pPamlNode)sibling)->dS;

		pTreeNode parent = ((pTreeNode)pnode)->parent;
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
		freenode((pTree)ptree,parent,NULL);	
	}
  
	Tree tmp;
	((pTreeNode)pnode)->parent = NULL;
	tmp.nlist = ((pTree)ptree)->nlist;
	tmp.root = (pTreeNode)pnode;
	tmp.count = ((pTree)ptree)->count;
	tree_traveral_prefix(&tmp, freenode );
}

int paml_tree_delete_node( pPamlTree ptree, pPamlNode pnode, tree_func_node freenode )
{
	if ( tree_is_root( (pTree)ptree, (pTreeNode)pnode) )
	{
		paml_tree_free(ptree); 
		return -1;
	}
  __paml_tree_delete_nodes(ptree,(pPamlNode)pnode,freenode);
  ((pPhylogenyNode)((pTree)ptree)->root)->branchlength = -1;
	tree_build_node_list((pTree)ptree);
  return 0;
}

/*
 * reg == 0 : delete
 * reg == 1 : remain
 */
int paml_tree_delete_nodes_by_reg( pPamlTree ptree, tree_func_node freenode )
{
	if ( ((pTree)ptree)->root->reg == 0 )
	{
		paml_tree_free(ptree); 
		return -1;
	}
	int i, bupdate = 0;
	pArrayList nlist = ((pTree)ptree)->nlist;
  
	for ( i = 0 ; i < nlist->size ; i++ )
	{
		pTreeNode pnode = (pTreeNode)nlist->array[i];
		if ( pnode == NULL ) continue;
		if ( ((pTree)ptree)->root == pnode ) {
      if (pnode->reg == 0) {
        // the only remaining node needs to be deleted, meaning the whole tree must be deleted
        if ( bupdate )
        {
          ((pPhylogenyNode)((pTree)ptree)->root)->branchlength = -1;
          tree_build_node_list((pTree)ptree);
        }
        paml_tree_free(ptree); 
        return -1;
      }
      else {
        continue;
      }
    }
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
				__paml_tree_delete_nodes(ptree,(pPamlNode)pnode,freenode);
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
		((pPhylogenyNode)((pTree)ptree)->root)->branchlength = -1;
		tree_build_node_list((pTree)ptree);
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

pTree paml_tree_split_tree(pTree ptree, pTreeNode pnode, tree_func_node freenode ) 
{
	if ( tree_is_root((pTree)ptree,(pTreeNode)pnode) ) return NULL;			
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
    ((pPamlNode)ptree->root)->N = 0;
    ((pPamlNode)ptree->root)->S = 0;
    ((pPamlNode)ptree->root)->dS = 0;
    ((pPamlNode)ptree->root)->dN = 0;
    ((pPamlNode)ptree->root)->omega = 0;
    ((pPamlNode)ptree->root)->accudS = 0;
	}
	else
	{
		((pPhylogenyNode)sibling)->branchlength += ((pPhylogenyNode)sibling->parent)->branchlength;
		((pPamlNode)sibling)->N += ((pPamlNode)sibling->parent)->N;
		((pPamlNode)sibling)->S += ((pPamlNode)sibling->parent)->S;
		((pPamlNode)sibling)->dN += ((pPamlNode)sibling->parent)->dN;
		((pPamlNode)sibling)->dS += ((pPamlNode)sibling->parent)->dS;
		((pPamlNode)sibling)->omega += ((pPamlNode)sibling->parent)->omega;
		((pPamlNode)sibling)->accudS = ((pPamlNode)sibling->parent)->accudS;
    
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
  ((pPamlTree)psubtree)->nei_gojobori_w = NULL;
  ((pPamlTree)psubtree)->nei_gojobori_dN = NULL;
  ((pPamlTree)psubtree)->nei_gojobori_dS = NULL;
  ((pPamlTree)psubtree)->gene_names = NULL;
  ((pPamlTree)psubtree)->omega = NULL;
  ((pPamlTree)psubtree)->prob = NULL;
  ((pPamlTree)psubtree)->nodelist = NULL;
  
	return psubtree;
}

