#include <cctype>
#include <algorithm>
#include <vector>
#include <sstream>
#include <cfloat>
#include <iomanip>
#include <fstream>

#include "error_model.h"
#include "cafe_commands.h"

extern "C" {
#include "cafe_shell.h"
	int __check_error_model_columnsums(pErrorStruct errormodel);
	void cafe_shell_free_errorstruct(pErrorStruct errormodel);
}

pErrorStruct cafe_shell_create_error_matrix_from_estimate(pErrorMeasure errormeasure);

struct to_lower {
	int operator() (int ch)
	{
		return std::tolower(ch);
	}
};

bool case_insensitive_equal(std::string s1, std::string s2)
{
	std::transform(s1.begin(), s1.end(), s1.begin(), to_lower());
	std::transform(s2.begin(), s2.end(), s2.begin(), to_lower());
	return s1 == s2;
}

std::vector<std::string> split(std::string str, char delimiter) 
{
	std::vector<std::string> result;
	std::stringstream ss(str); // Turn the string into a stream.
	std::string tok;

	while (getline(ss, tok, delimiter)) {
		result.push_back(tok);
	}

	return result;
}

pErrorStruct get_error_model(pCafeFamily family, std::string filename)
{
	pErrorStruct result = NULL;
	if (family->errors) {
		for (int i = 0; i<family->errors->size; i++) {
			pErrorStruct error = (pErrorStruct)family->errors->array[i];
			if (case_insensitive_equal(error->errorfilename, filename)) {
				result = error;
				break;
			}
		}
	}

	return result;
}

std::ostream& operator<<(std::ostream& ost, ErrorStruct& errormodel)
{
	int i, j;

	ost << "maxcnt:" << errormodel.maxfamilysize << "\n";
	ost << "cntdiff";
	for (j = errormodel.fromdiff; j <= errormodel.todiff; j++) {
		ost << " " << j;
	}
	ost << "\n";

	ost << std::setw(2) << std::setprecision(2);
	for (j = 0; j <= errormodel.maxfamilysize; j++) {
		ost << j;
		for (i = errormodel.fromdiff; i <= errormodel.todiff; i++) {
			if (0 <= i + j && i + j <= errormodel.maxfamilysize) {
				ost << " " << errormodel.errormatrix[i + j][j];	// conditional probability of measuring i+j when true count is j        
			}
			else {
				ost << " #nan";
			}
		}
		ost << "\n";
	}

	return ost;
}

std::istream& operator>>(std::istream& ifst, ErrorStruct& errormodel)
{
	std::string line;
	if (!std::getline(ifst, line))
	{
		throw std::runtime_error("Empty file");
	}
	std::vector<std::string> data = split(line, ' ');
	std::vector<std::string> max = split(data.at(0), ':');
	errormodel.maxfamilysize = atoi(max[1].c_str());

	if (std::getline(ifst, line)) {
		std::vector<std::string> data = split(line, ' ');
		errormodel.fromdiff = atoi(data.at(1).c_str());
		errormodel.todiff = atoi(data.at(data.size() - 1).c_str());
	}
	errormodel.errormatrix = (double**)memory_new_2dim(errormodel.maxfamilysize + 1, errormodel.maxfamilysize + 1, sizeof(double));

	int i = 0;
	int j = 0;
	while (std::getline(ifst, line))
	{
		data = split(line, ' ');
		int col1 = atoi(data[0].c_str());
		if ((int)data.size() == (errormodel.todiff - errormodel.fromdiff) + 2) {
			while (j && j < col1) {
				// copy previous line's error model for missing lines. 
				for (i = errormodel.fromdiff; i <= errormodel.todiff; i++) {

					if (i + j >= 0 && i + j <= errormodel.maxfamilysize) {
						errormodel.errormatrix[i + j][j] = errormodel.errormatrix[i + j - 1][j - 1];
					}
				}
				i++;
			}
			// read error model and save in matrix row
			int k = 1;  // k is file column index
			for (i = errormodel.fromdiff; i <= errormodel.todiff; i++) {
				assert(j == col1);
				if (i + j >= 0 && i + j <= errormodel.maxfamilysize) {
					errormodel.errormatrix[i + j][j] = atof(data[k].c_str());  // conditional probability of measuring i+j when true count is j
				}
				k++;
			}
			j++;
		}
	}
	while (j && j <= errormodel.maxfamilysize) {
		// copy previous line's error model for missing lines till the end of matrix. 
		for (i = errormodel.fromdiff; i <= errormodel.todiff; i++) {
			if (i + j >= 0 && i + j <= errormodel.maxfamilysize) {
				errormodel.errormatrix[i + j][j] = errormodel.errormatrix[i + j - 1][j - 1]; // conditional probability of measuring i+j when true count is j
			}
		}
		j++;
	}

	return ifst;
}

void init_error_ptr(pCafeFamily family, pCafeTree pTree, pErrorStruct errormodel, std::string speciesname)
{
	if (family->error_ptr == NULL) {
		family->error_ptr = (pErrorStruct *)memory_new(family->num_species, sizeof(pErrorStruct));
	}
	if (!speciesname.empty()) {
		for (int i = 0; i<family->num_species; i++) {
			if (case_insensitive_equal(family->species[i], speciesname)) {
				family->error_ptr[i] = errormodel;
				pCafeNode pcnode = (pCafeNode)pTree->super.nlist->array[family->index[i]];
				pcnode->errormodel = errormodel;
				break;
			}
		}
	}
	else { // '-all' specified instead of speciesname 
		for (int i = 0; i<family->num_species; i++) {
			family->error_ptr[i] = errormodel;
			pCafeNode pcnode = (pCafeNode)pTree->super.nlist->array[family->index[i]];
			pcnode->errormodel = errormodel;
		}
	}

}

int set_error_matrix_from_file(pCafeFamily family, pCafeTree pTree, family_size_range& range, std::string filename, std::string speciesname)
{
	// check if error model for filename already exists 
	pErrorStruct errormodel = get_error_model(family, filename);

	if (errormodel == NULL)
	{
		// allocate new errormodel
		errormodel = (pErrorStruct)calloc(1, sizeof(ErrorStruct));
		errormodel->errorfilename = strdup(filename.c_str());
		std::ifstream ifst(filename.c_str());
		if (!ifst)
			throw io_error("errormodel", filename, false);

		ifst >> *errormodel;

		if (errormodel->maxfamilysize < range.max) {
			errormodel->maxfamilysize = range.max;
		}
		// now make sure that columns of the error matrix sums to one.
		__check_error_model_columnsums(errormodel);

		if (family->errors == NULL) {
			family->errors = arraylist_new(family->num_species);
		}
		arraylist_add(family->errors, errormodel);

	}
	init_error_ptr(family, pTree, errormodel, speciesname);
	return 0;
}

int remove_error_model(pCafeFamily family, pCafeTree pcafe, std::string species_name)
{
	int i = 0;
	if (family->errors) {
		assert(family->error_ptr != NULL);	// errors and error_ptr move in lockstep

		for (i = 0; i<family->num_species; i++) {
			if (case_insensitive_equal(family->species[i], species_name)) {
				family->error_ptr[i] = NULL;
				pCafeNode pcnode = (pCafeNode)pcafe->super.nlist->array[family->index[i]];
				pcnode->errormodel = NULL;
				break;
			}
		}
	}
	return 0;
}

void free_error_model(pCafeFamily family, pCafeTree pcafe)
{
	int i;
	for (i = 0; i<family->num_species; i++) {
		remove_error_model(family, pcafe, family->species[i]);
	}
	if (family->errors) {
		arraylist_free(family->errors, (freefunc)cafe_shell_free_errorstruct);
		family->errors = NULL;
	}
	if (family->error_ptr) {
		memory_free(family->error_ptr);
		family->error_ptr = NULL;
	}
}


