 #include<memalloc.h>
#include<assert.h>

void* memory_new(size_t count, size_t size )
{
  void* data = calloc(count,size);
	if ( data == NULL )
	{
		fprintf(stderr,"(%s) Error in memory allocation: %d: %s\n", __FUNCTION__, errno, strerror(errno) );
		exit(errno);
	}
	return data;  
}

void** memory_new_2dim(int row, int col, int size )
{
	void** data =(void**)calloc(row,sizeof(void*));			
	if ( data == NULL )
	{
		fprintf(stderr,"(%s) Error in memory allocation: %d: %s\n", __FUNCTION__, errno, strerror(errno) );
		perror(__FUNCTION__);
		exit(errno);
	}
	int i;
	for ( i = 0 ; i < row; i++  )
	{
		data[i] = calloc(col,size);
		if ( data[i] == NULL )
		{
			fprintf(stderr,"(%s) Error in memory allocation: %d: %s\n", __FUNCTION__, errno, strerror(errno) );
			perror(__FUNCTION__);
			exit(errno);
		}
	}
	/*int r;
	for ( r = 0 ; r < row ; r++ )
	{
		memset(data[r], 0, col * size );
	}*/
	return data;
}


void* memory_new_with_init(int count, int size, void* init)
{
	void* data = memory_new(count,size);
	memcpy(data,init,count*size);
	return data;
}

void* memory_realloc_with_data(void** data, int old, int newsize, int size, func_memory_free func, int element_size )
{
	int i;
	if ( old == newsize ) return data;
	if ( old > newsize && func )
	{
		for ( i = newsize ; i < old ; i++ )
		{
			func(data[i]);	
		}
	}
	data = (void**)realloc(data, newsize * size );
	if ( newsize > old && element_size > 0 )
	{
		for ( i = old ; i < newsize ; i++ )
		{
			data[i] = memory_new(1,element_size);
		}
	}
	return data;
}

void* memory_realloc(void* data, int count, int size )
{
	if ( (data = realloc(data, count*size)) == NULL )
	{
		fprintf(stderr,"Error in memory reallocation with clear: %d\n", errno );
		exit(errno);
	}
	return data;
}

void memory_free(void* data)
{
	assert(data != NULL);
	free(data);
}

void memory_free_2dim(void** data, int row, int col, func_memory_free func )
{
	int r;
	assert(data != NULL);
	if ( func )
	{
		for ( r = 0 ; r < row ; r++ )
		{
			func(data[r]);
            data[r] = NULL;
		}
	}
	else
	{
		for ( r = 0 ; r < row ; r++ )
		{
			free(data[r]);
            data[r] = NULL;
		}
	}
	free(data);
}

void** memory_copy_2dim(void** dst, void** src, int row, int col, int elem_size )
{
	int r; 
	for ( r = 0 ; r < row ; r++ )
	{
		memcpy(dst[r], src[r], col*elem_size );
	}
	return dst;
}

