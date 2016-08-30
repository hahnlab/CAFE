#ifndef __REG_EXPRESS_H__
#define __REG_EXPRESS_H__

#include<utils.h>
#include<utils_string.h>
#include<regex.h>

extern pArrayList regex_split(char* pattern, char* string);
extern size_t regex_match(char* pattern, char* string, int opt, regmatch_t* match);
extern pArrayList regex_match_group(char* pattern, char* string, int opt, regmatch_t* match);

#endif 
