#include<stdio.h>
#include<memalloc.h>
#include<utils_string.h>

typedef struct
{
	char* name;
	char* value;
}CGIParam;

typedef CGIParam* pCGIParam; 

pCGIParam cgi_param;
int num_cgi_params;

void cgi_parse_parameter_of_get()
{
	cgi_param = NULL;
	char* str = getenv("QUERY_STRING");
	if ( str == NULL  || strlen(str) == 0 ) return; 
	pArrayList ap = string_pchar_split( str, '&' );
	cgi_param = (pCGIParam) memory_new(ap->size,sizeof(CGIParam));
	num_cgi_params = ap->size;
	int i;
	for ( i = 0 ; i < ap->size ; i++ )
	{
		cgi_param[i].name = (char*)ap->array[i];
		char* eq = index((char*)ap->array[i],'=');
		if ( eq ) 
		{
			*eq = '\0';
			cgi_param[i].value = (char*)(eq+1);
		}
	}
}

void cgi_free()
{
	int i;
	for ( i = 0 ; i < num_cgi_params; i++ )
	{
		memory_free( cgi_param[i].name );
		cgi_param[i].name = NULL;
	}
	memory_free( cgi_param );
	cgi_param = NULL;
}

void cgi_set_type(char* type)
{
	printf("Content-type: %s\n\n", type );
}

void cgi_start_html(char* title)
{
	printf("<html>\n<head>\t<title>%s</title>\n</head>\n", title );
}

void cgi_end_html()
{
	printf("</html>\n");
}

char* cgi_get_parameter(char* name)
{
	int i;
	for ( i = 0 ; i < num_cgi_params ; i++ )
	{
		if ( strcmp( cgi_param[i].name, name ) == 0 ) return cgi_param[i].value;
	}
	return NULL;
}
