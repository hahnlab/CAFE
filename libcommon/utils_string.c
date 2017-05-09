
#include <stdlib.h>
#include "utils_string.h"
#include "memalloc.h"

pString string_new_step(size_t step)
{
	pString pstr = (pString)memory_new(1,sizeof(String));
	pstr->buf = (char*)memory_new(step+1,sizeof(char));
	pstr->buf[0] = '\0';
	pstr->length = 0;
	pstr->alloc_size = step;
	return pstr;
}

pString string_new()
{
	return string_new_step(STRING_STEP_SIZE);
}

pString string_new_with_string(const char* newstr )
{
	size_t len = strlen(newstr);
	pString pstr = string_new_step(len+1);		
	strcpy(pstr->buf,newstr);
	pstr->alloc_size = len;
	pstr->length = len;
	return pstr;
}

void string_free(pString pstr)
{
	memory_free(pstr->buf);
	pstr->buf = NULL;
	memory_free(pstr);
	pstr = NULL;
}

void string_free_without_data(pString pstr)
{
	memory_free(pstr);
	pstr = NULL;
}

void string_reset(pString pstr)
{
	pstr->length = 0;
	pstr->buf[0] = '\0';
}

void string_add(pString pstr, const char* add)
{
	size_t len = strlen(add);
	if ( pstr->alloc_size < len + pstr->length )	
	{
		pstr->alloc_size += len + STRING_STEP_SIZE;
		pstr->buf = (char*)realloc(pstr->buf, pstr->alloc_size + 1 );
		if ( pstr->buf == NULL )
		{
			print_error(__FILE__,(char*)__FUNCTION__,__LINE__,
					    "memory reallocation of String->buf");		
		}
	}
	strcat(pstr->buf,add);
	pstr->length += len;
}

void string_fadd(pString pstr, const char* msg, ... )
{
	char buf[STRING_BUF_SIZE];
	va_list ap;
	va_start(ap,msg);
	vsprintf(buf,msg,ap);
	va_end(ap);

	string_add(pstr,buf);
}

void string_trim(pString pstr)
{
	pstr->buf = (char*)realloc(pstr->buf, pstr->length + 1 );
	if ( pstr->buf == NULL )
	{
		print_error(__FILE__,(char*)__FUNCTION__,__LINE__,
					 "memory reallocation of String->buf");		
	}
	pstr->alloc_size = pstr->length;
}

void string_chomp(pString pstr)
{
	string_pchar_chomp(pstr->buf);	
	size_t newlen = strlen(pstr->buf);
	pstr->length = newlen;
}

char* string_get(pString pstr)
{
	return pstr->buf;
}

int is_letter(char c)
{
	return ( 'A' <= c && 'Z' >= c ) || ( 'a' <= c && 'z' >= c );
}

int is_space(char c)
{
	return c == ' ' || c == '\t'  || c == '\n'  || c == '\r';
}

char* string_pchar_chomp(char* pstr)
{
	long from, end;
	for ( from = 0 ; pstr[from] ; from++ )
	{
		if ( !is_space(pstr[from]) ) break;;
	}
	for ( end = strlen(pstr) - 1 ; end >= 0 ; end-- )
	{
		if ( !is_space(pstr[end]) ) break;;
	}
	size_t i, j;
    if (from <= end) {
        if (from > 0) {
            for ( i = from, j = 0 ; i <= end ; i++, j++ )
            {
                pstr[j] = pstr[i];
            }
        }
        pstr[end-from+1] = '\0';
    }
	return pstr;
}

int string_pchar_cmp_ignore_case(char* cmp1, char* cmp2)
{
	int i;
	size_t len = strlen(cmp2);
    size_t len1 = strlen(cmp1);
    if (len != len1) {
        return 0;
    }
	for ( i = 0 ; i < len ; i++ )
	{
		if ( is_letter(cmp1[i]) && is_letter(cmp2[i]) )
		{
			if ( (cmp1[i] & 0x5F) != (cmp2[i] & 0x5F) ) return 0;
		}
		else
		{
			if ( cmp1[i] != cmp2[i]) return 0;
		}
	}
	return 1;
}

pArrayList string_pchar_split(char* buf, char delim)
{
	pArrayList psplit = arraylist_new(20);
	char* begin = buf;
	char* str;
	size_t len;
	while( *buf )
	{
		if ( *buf == delim )
		{
			*buf++ = '\0';
			len = strlen(begin);
			str = NULL;
			if ( len )
			{
				str = (char*)memory_new( len+1, sizeof(char) );
				strcpy(str,begin);
			}
			arraylist_add(psplit,str);
			begin = buf;
		}
		else
		{
			buf++;
		}
	}
	len = strlen(begin);
	str = NULL;
	if ( len )
	{
		str = (char*)memory_new( len+1, sizeof(char) );
		strcpy(str,begin);
	}
	arraylist_add(psplit,str);
	return psplit;
}

pArrayList string_pchar_space_split(char* buf)
{
	pArrayList parg = arraylist_new(10);

    int i = 0;
	while( *buf == ' ' || *buf == '\t' ) buf++;
	char* start = buf;
	while( *buf )
	{
		if ( *buf == ' ' || *buf == '\t' || *buf == '\n' )
		{
			*buf++ = '\0';
			if ( i > 0 )
			{
				arraylist_add( parg, start );
			}
			start = buf;
			i = 0;
		}
		else
		{
			i++;
			buf++;
		}
	}
	if ( i > 0 )
	{
		arraylist_add( parg, start );
	}
	return parg;
}

void string_pchar_join_double(char* rtn, const char* sp, int argc, double* values)
{
	int i;
	char buf[STRING_STEP_SIZE];
	for ( i = 0 ; i < argc ; i++ )
	{
		sprintf(buf,"%15.14lf", values[i] );
		strcat(rtn,buf);
		if ( i < argc - 1 ) strcat(rtn,sp);
	}
}

void string_pchar_join(char* buf, char* stuff, int num, char** list)
{
	int i;
	strcpy(buf,list[0]);
	for ( i = 1 ; i < num ; i++ )
	{
		if ( stuff ) strcat(buf, stuff);
		strcat(buf, list[i]);	
	}
}

pString string_join(const char* stuff, int num, char** list)
{
	int i;
	pString pstr = string_new_with_string(list[0]);
	for( i = 1 ; i < num ; i++ )
	{
		string_add(pstr,stuff);
		string_add(pstr,list[i]);
	}
	return pstr;
}
