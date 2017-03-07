#include "regexpress.h"
#include "memalloc.h"

pArrayList regex_split(char* pattern, char* string)
{
	int err;
	pArrayList pal = arraylist_new(-1);
	regex_t re;
	regmatch_t match;

	if ( (err=regcomp(&re,pattern, REG_EXTENDED )) )
	{
		char buf[STRING_BUF_SIZE];
		regerror( err, &re, buf, STRING_BUF_SIZE );
		fprintf( stderr, "%s", buf );
		return NULL; 
	}

	while( string && !(err=regexec(&re,string,1,&match,0)) )
	{
		size_t len = (size_t)match.rm_so;
		char* str = (char*)memory_new(len+1,sizeof(char));
		strncpy(str,string,len);
		str[len] = '\0';
		arraylist_add(pal,str);
		string += match.rm_eo;
	}
	
	if ( string )
	{
		size_t len = strlen(string);
		char* str = (char*)memory_new(len+1,sizeof(char));
		strcpy(str,string);
		arraylist_add(pal,str);
	}

	if ( err && err != REG_NOMATCH )
	{
		char buf[STRING_BUF_SIZE];
		regerror( err, &re, buf, STRING_BUF_SIZE );
		fprintf( stderr, "%s", buf );
		return NULL; 
	}
	regfree(&re);
	return pal;
}

size_t regex_match(char* pattern, char* string, int opt, regmatch_t* match)
{
	int err;
	regex_t re;

	if ( (err=regcomp(&re,pattern, REG_EXTENDED | opt )) )
	{
		char buf[STRING_BUF_SIZE];
		regerror( err, &re, buf, STRING_BUF_SIZE );
		fprintf( stderr, "%s", buf );
		return 0;
	}
	err=regexec(&re,string,1,match,0);
	regfree(&re);
	return err ? -err : (size_t)match->rm_eo;
}

pArrayList regex_match_group(char* pattern, char* string, int opt, regmatch_t* match)
{
	int grp_cnt = 0;
	int grp_stack = 0;
	size_t pat_len = strlen(pattern);
	char* tmp_pattern = (char*) memory_new( pat_len, sizeof(char) );
	strcpy(tmp_pattern, pattern);
		
	int i, j;
	for ( i = j = 0 ; i < pat_len ; i++ )
	{
		if ( i != 0 && pattern[i-1] != '\\' && pattern[i] == '(' ) 
		{
			grp_cnt++;
			grp_stack++;
		}
		else if ( i != 0 && pattern[i-1] != '\\' && pattern[i] == ')' ) 
		{
			grp_stack--;
		}
		else
		{
			tmp_pattern[j++] = pattern[i];
		}	
	}
	if( grp_stack != 0 )
	{
		fprintf( stderr, "Mismatch parenthesis in %s\n", pattern );
		memory_free(tmp_pattern);
		tmp_pattern = NULL;
		return NULL;	
	}

	memory_free(tmp_pattern);
	tmp_pattern = NULL;
	return NULL;
}
