#ifndef __UTILS_STRING_H__
#define __UTILS_STRING_H__

#include<utils.h>
#include<strings.h>
#include<stdio.h>
#include<memalloc.h>

#define STRING_STEP_SIZE 2048
#define STRING_BUF_SIZE 131072

typedef struct
{
	char* buf;
	size_t length;
	size_t alloc_size;
}String;
typedef String* pString;

extern pString string_new_step(size_t step);
extern pString string_new();
extern pString string_new_with_string(const char* newstr );
extern void string_free(pString pstr);
extern void string_free_without_data(pString pstr);
extern void string_reset(pString pstr);
extern void string_add(pString pstr, const char* add);
extern void string_fadd(pString pstr, const char* msg, ...);
extern void string_trim(pString pstr);
extern char* string_get(pString pstr);

extern char* string_pchar_chomp(char* pstr);
extern int string_pchar_cmp_ignore_case(char* cmp1, char* cmp2);
extern pArrayList string_pchar_split(char* buf, char delim);
extern pArrayList string_pchar_space_split(char* buf);
extern void string_pchar_join_double(char* rtn, char* sp, int argc, double* values);
extern void string_pchar_join(char* buf, char* stuff, int num, char** list);
extern pString string_join(const char* stuff, int num, char** list);

#endif
