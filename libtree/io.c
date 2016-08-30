#include"io.h"
#include<utils_string.h>
#include<dirent.h>

char* file_read_all(char* fname)
{
	struct stat fstat;
	stat(fname,&fstat);
	size_t size = (size_t)fstat.st_size;
	char* buf = (char*) calloc(size+1,sizeof(char));
	if ( buf == NULL )
	{
		print_error(__FILE__,(char*)__FUNCTION__,__LINE__,
		             "memory allocation of buf");		
	}

	int fh = open( fname, O_RDONLY ); 	
	if ( fh == -1 )
	{
		print_error(__FILE__,(char*)__FUNCTION__,__LINE__,
		             "Cannot open file: %s", fname);		
	}
	ssize_t i;
	if ( (i=read(fh,buf,size)) <= 0 )
	{
		if ( i == 0 )
		{
			print_error(__FILE__,(char*)__FUNCTION__,__LINE__,
		             "Empty file, %s", fname);		
		}
		else
		{
			print_error(__FILE__,(char*)__FUNCTION__,__LINE__,
		             "Read error of file, %s", fname);		
		}
	}
	close(fh);

	buf[size] = '\0';
	return buf;
}

size_t file_read(pString pstr, FILE* fp)
{
	char buf[STRING_BUF_SIZE];
	size_t s;
	string_reset(pstr);
	while( (s=fread( buf, 1, STRING_BUF_SIZE-1, fp )) )
	{
		buf[s] = '\0';	
		string_add(pstr, buf);
	}
	return pstr->length;
}

size_t file_read_line(pString pstr, FILE* fp)
{
	char buf[STRING_BUF_SIZE];
	string_reset(pstr);
	while( fgets(buf,STRING_BUF_SIZE,fp))
	{
		string_add(pstr, buf);
		size_t len = strlen(buf);
		if ( buf[len-1] == '\n' ) break;
	}
	return pstr->length;
}

pArrayList dir_file_list(char* dir, char* ext)
{
	DIR* d = opendir(dir);
	struct dirent* dp;
	if ( d == NULL ) return NULL;

	pArrayList pal = arraylist_new(-1);	
	while( (dp=readdir(d)) != NULL )
	{
		if ( ext )
		{
			char* point = rindex(dp->d_name, '.' );
			if ( point == NULL ) continue;
			point++;
			if ( strcmp( point, ext ) ) continue;
		}
		char* name = (char*) memory_new( strlen(dp->d_name)+1,sizeof(char));
		strcpy(name, dp->d_name);
		arraylist_add(pal,name);		
	}
	closedir(d);
	return pal;
}


