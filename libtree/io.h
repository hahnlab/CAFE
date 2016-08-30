#ifndef __IO_H__
#define __IO_H__

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<sys/uio.h>
#include<fcntl.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<unistd.h>
#include<utils.h>
#include<utils_string.h>

extern char* file_read_all(char* fname);
extern size_t file_read_line(pString pstr, FILE* fp);
extern size_t file_read(pString pstr, FILE* fp);
extern pArrayList dir_file_list(char* dir, char* ext);

#endif
