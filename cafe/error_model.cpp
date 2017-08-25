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
	int cafe_shell_read_error_double_measure(const char* error1, const char* error2, int** observed_pairs, int maxFamilySize);
	double __loglikelihood_pairs_from_double_measure(double* parameters, void* args);
	int cafe_shell_read_error_true_measure(const char* errorfile, const char* truefile, int** observed_pairs, int maxFamilySize);
	double __loglikelihood_pairs_from_true_measure(double* parameters, void* args);
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
  int file_row_count = atoi(max[1].c_str());
	errormodel.maxfamilysize = std::max(errormodel.maxfamilysize, file_row_count);

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
    errormodel->maxfamilysize = range.max;
		std::ifstream ifst(filename.c_str());
		if (!ifst)
			throw io_error("errormodel", filename, false);

		ifst >> *errormodel;

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

int cafe_shell_read_freq_from_measures(const char* file1, const char* file2, int* sizeFreq, int maxFamilySize)
{
	int i = 0;
	char buf1[STRING_BUF_SIZE];
	char buf2[STRING_BUF_SIZE];

	FILE* fpfile1 = fopen(file1, "r");
	if (fpfile1 == NULL)
	{
		fprintf(stderr, "Cannot open file: %s\n", file1);
		return -1;
	}
	if (fgets(buf1, STRING_BUF_SIZE, fpfile1) == NULL)
	{
		fclose(fpfile1);
		fprintf(stderr, "Empty file: %s\n", file1);
		return -1;
	}
	FILE* fpfile2 = NULL;
	if (file2) {
		fpfile2 = fopen(file2, "r");
		if (fpfile2 == NULL)
		{
			fprintf(stderr, "Cannot open file: %s\n", file2);
			return -1;
		}
		if (fgets(buf2, STRING_BUF_SIZE, fpfile2) == NULL)
		{
			fclose(fpfile2);
			fprintf(stderr, "Empty file: %s\n", file2);
			return -1;
		}
	}

	string_pchar_chomp(buf1);
	pArrayList data1 = string_pchar_split(buf1, '\t');
	arraylist_free(data1, NULL);

	// count frequency of family sizes
	int line1 = 0;
	// int maxFamilySize = cafe_param->family_size.max;
	int data1colnum = 0;
	while (fgets(buf1, STRING_BUF_SIZE, fpfile1))
	{
		string_pchar_chomp(buf1);
		pArrayList data1 = string_pchar_split(buf1, '\t');
		for (i = 2; i<data1->size; i++) {
			int size1 = atoi((char*)data1->array[i]);
			sizeFreq[size1]++;
			if (size1 > maxFamilySize) {
				maxFamilySize = size1;
			}
		}
		data1colnum = data1->size;
		arraylist_free(data1, NULL);
		line1++;
	}
	if (fpfile2)
	{
		int line2 = 0;
		string_pchar_chomp(buf2);
		pArrayList data2 = string_pchar_split(buf2, '\t');
		if (data1colnum != data2->size) {
			fprintf(stderr, "file: the number of columns do not match between the two files\n");
			return -1;
		}
		arraylist_free(data2, NULL);

		while (fgets(buf2, STRING_BUF_SIZE, fpfile2))
		{
			string_pchar_chomp(buf2);
			pArrayList data2 = string_pchar_split(buf2, '\t');
			for (i = 2; i<data2->size; i++) {
				int size2 = atoi((char*)data2->array[i]);
				sizeFreq[size2]++;
				if (size2 > maxFamilySize) {
					maxFamilySize = size2;
				}
			}
			arraylist_free(data2, NULL);
			line2++;
		}
		if (line1 != line2) {
			fprintf(stderr, "ERROR: the number of lines do not match between the two files\n");
			return -1;
		}
	}
	return maxFamilySize;
}




int cafe_shell_read_error_double_measure(const char* error1, const char* error2, int** observed_pairs, int maxFamilySize)
{
	int i = 0;
	int j = 0;
	char buf1[STRING_BUF_SIZE];
	char buf2[STRING_BUF_SIZE];

	FILE* fperror1 = fopen(error1, "r");
	if (fperror1 == NULL)
	{
		fprintf(stderr, "Cannot open file: %s\n", error1);
		return -1;
	}
	if (fgets(buf1, STRING_BUF_SIZE, fperror1) == NULL)
	{
		fclose(fperror1);
		fprintf(stderr, "Empty file: %s\n", error1);
		return -1;
	}
	FILE* fperror2 = fopen(error2, "r");
	if (fperror2 == NULL)
	{
		fprintf(stderr, "Cannot open file: %s\n", error2);
		return -1;
	}
	if (fgets(buf2, STRING_BUF_SIZE, fperror2) == NULL)
	{
		fclose(fperror2);
		fprintf(stderr, "Empty file: %s\n", error2);
		return -1;
	}


	// now compare two files and count pairs.
	while (fgets(buf1, STRING_BUF_SIZE, fperror1))
	{
		if (fgets(buf2, STRING_BUF_SIZE, fperror2) != NULL) {
			string_pchar_chomp(buf1);
			pArrayList data1 = string_pchar_split(buf1, '\t');
			string_pchar_chomp(buf2);
			pArrayList data2 = string_pchar_split(buf2, '\t');
			if (strcmp((char*)data1->array[1], (char*)data2->array[1]) != 0) {
				fprintf(stderr, "ERROR: the family IDs in each line do not match between the two files\n");
				return -1;
			}
			// check pairs
			for (i = 2; i<data1->size; i++) {
				int size1 = atoi((char*)data1->array[i]);
				int size2 = atoi((char*)data2->array[i]);
				observed_pairs[size1][size2]++;
			}
			arraylist_free(data1, NULL);
			arraylist_free(data2, NULL);
		}
	}

	// now make triangle matrix by merging i,j and j,i
	for (i = 0; i <= maxFamilySize; i++) {
		for (j = 0; j<i; j++) {
			observed_pairs[j][i] += observed_pairs[i][j];
			observed_pairs[i][j] = 0;
		}
	}


	return 0;
}


pErrorMeasure estimate_error_double_measure(std::ostream& log, const char* error1, const char* error2, int b_symmetric, int max_diff, int b_peakzero, int max_FamilySize)
{
	int i;
	int* sizeFreq = (int *)memory_new(10000, sizeof(int));
	int maxFamilySize = cafe_shell_read_freq_from_measures(error1, error2, sizeFreq, max_FamilySize);
	if (maxFamilySize < 0) {
		fprintf(stderr, "ERROR: failed to read freqeuncy from measurement files\n");
	}
	// get size probability distribution
	int sizeTotal = 0;
	for (i = 0; i <= maxFamilySize; i++) {
		sizeTotal += sizeFreq[i] + 1;
		if (sizeTotal < 0) {
			fprintf(stderr, "ERROR: total freqeuncy is less than zero\n");
		}
	}
	double* sizeDist = (double*)memory_new(maxFamilySize + 1, sizeof(double));
	for (i = 0; i <= maxFamilySize; i++) {
		sizeDist[i] = (sizeFreq[i] + 1) / (double)sizeTotal;
		if (sizeDist[i] < 0) {
			fprintf(stderr, "ERROR: freqeuncy is less than zero\n");
		}
	}


	int** observed_pairs = (int**)memory_new_2dim(maxFamilySize + 1, maxFamilySize + 1, sizeof(int)); // need space for zero
	int retval = cafe_shell_read_error_double_measure(error1, error2, observed_pairs, maxFamilySize);
	if (retval < 0) {
		fprintf(stderr, "ERROR: failed to count pairs from measurement files\n");
	}

	// set up parameters for ML
	pErrorMeasure error = (pErrorMeasure)memory_new(1, sizeof(ErrorMeasure));
	error->sizeDist = sizeDist;
	error->maxFamilySize = maxFamilySize;
	error->pairs = observed_pairs;
	error->b_symmetric = b_symmetric;
	error->b_peakzero = b_peakzero;
	if (b_symmetric) {
		// symmetric model (diff == number)
		error->model_parameter_diff = max_diff;
		error->model_parameter_number = max_diff + 1;
	}
	else {
		// asymmetric model (diff*2 == number)
		error->model_parameter_diff = max_diff;
		error->model_parameter_number = 2 * max_diff + 1;
	}


	// now estimate the misclassification rate 
	int max_runs = 100;
	int converged = 0;
	int runs = 0;
	double minscore = DBL_MAX;
	double* parameters = (double *)memory_new(error->model_parameter_number, sizeof(double));
	double* bestrun_parameters = (double *)memory_new(error->model_parameter_number, sizeof(double));

	do {
		pFMinSearch pfm;
		double* sorted_params = (double *)memory_new_with_init(error->model_parameter_number, sizeof(double), parameters);
		for (i = 0; i<error->model_parameter_number; i++) {
			sorted_params[i] = unifrnd() / (double)error->model_parameter_number;
		}
		qsort(sorted_params, error->model_parameter_number, sizeof(double), comp_double);
		if (error->b_symmetric) {
			int j = 0;
			for (i = error->model_parameter_number - 1; i >= 0; i--) {
				parameters[j++] = sorted_params[i];
			}
		}
		else {
			int j = error->model_parameter_number - 1;
			parameters[error->model_parameter_diff] = sorted_params[j--];
			for (i = 1; i <= error->model_parameter_diff; i++) {
				parameters[error->model_parameter_diff - i] = sorted_params[j--];
				parameters[error->model_parameter_diff + i] = sorted_params[j--];
			}
		}
		pfm = fminsearch_new_with_eq(__loglikelihood_pairs_from_double_measure, error->model_parameter_number, error);
		pfm->tolx = 1e-9;
		pfm->tolf = 1e-9;
		fminsearch_min(pfm, parameters);
		double *re = fminsearch_get_minX(pfm);
		for (i = 0; i < error->model_parameter_number; i++) parameters[i] = re[i];
		log << std::endl << "Misclassification Matrix Search Result: (" << pfm->iters << " iterations)\n";
		log << "Score: " << *pfm->fv << std::endl;

		if (runs > 0) {
			if (!isnan(*pfm->fv) && !isinf(*pfm->fv) && abs(minscore - (*pfm->fv)) < pfm->tolf) {
				converged = 1;
			}
		}
		if (pfm->iters < pfm->maxiters) {
			if (*pfm->fv < minscore) {
				minscore = *pfm->fv;
				memcpy(bestrun_parameters, parameters, (error->model_parameter_number) * sizeof(double));
			}
			runs++;
		}
		/*        else {
		cafe_log(param,"what went wrong?\n");
		fminsearch_min(pfm, parameters);
		}*/
		fminsearch_free(pfm);
	} while (!converged && runs<max_runs);

	if (converged) {
		log << "score converged in " << runs << " runs." << std::endl;
	}
	else {
		log << "score failed to converge in " << max_runs << " runs." << std::endl;
		log << "best score: " << minscore << std::endl;
	}
	memory_free(parameters);
	error->estimates = bestrun_parameters;

	//memory_free(error);           // we are going to return these values
	memory_free_2dim((void**)observed_pairs, maxFamilySize + 1, maxFamilySize + 1, NULL);
	memory_free(sizeFreq);
	return error;
}

pErrorMeasure estimate_error_true_measure(std::ostream& log, const char* errorfile, const char* truefile, int b_symmetric, int max_diff, int b_peakzero, int max_family_size)
{
	int i;

	int* sizeFreq = (int *)memory_new(10000, sizeof(int));
	int maxFamilySize = cafe_shell_read_freq_from_measures(truefile, errorfile, sizeFreq, max_family_size);
	if (maxFamilySize < 0) {
		fprintf(stderr, "ERROR: failed to read freqeuncy from measurement files\n");
	}
	// get size probability distribution
	int sizeTotal = 0;
	for (i = 0; i <= maxFamilySize; i++) {
		sizeTotal += sizeFreq[i] + 1;
	}
	double* sizeDist = (double*)memory_new(maxFamilySize + 1, sizeof(double));
	for (i = 0; i <= maxFamilySize; i++) {
		sizeDist[i] = (sizeFreq[i] + 1) / (double)sizeTotal;
	}


	int** observed_pairs = (int**)memory_new_2dim(maxFamilySize + 1, maxFamilySize + 1, sizeof(int)); // need space for zero
	int retval = cafe_shell_read_error_true_measure(errorfile, truefile, observed_pairs, maxFamilySize);
	if (retval < 0) {
		fprintf(stderr, "ERROR: failed to count pairs from measurement files\n");
	}

	// set up parameters for ML
	pErrorMeasure error = (pErrorMeasure)memory_new(1, sizeof(ErrorMeasure));
	error->sizeDist = sizeDist;
	error->maxFamilySize = maxFamilySize;
	error->pairs = observed_pairs;
	error->b_symmetric = b_symmetric;
	error->b_peakzero = b_peakzero;
	if (b_symmetric) {
		// symmetric model (diff == number)
		error->model_parameter_diff = max_diff;
		error->model_parameter_number = max_diff + 1;
	}
	else {
		// asymmetric model (diff*2 == number)
		error->model_parameter_diff = max_diff;
		error->model_parameter_number = 2 * max_diff + 1;
	}

	// now estimate the misclassification rate 
	int max_runs = 100;
	int converged = 0;
	int runs = 0;
	double minscore = DBL_MAX;
	double* parameters = (double *)memory_new(error->model_parameter_number, sizeof(double));
	double* bestrun_parameters = (double *)memory_new(error->model_parameter_number, sizeof(double));

	do {
		pFMinSearch pfm;
		double* sorted_params = (double *)memory_new_with_init(error->model_parameter_number, sizeof(double), parameters);
		for (i = 0; i<error->model_parameter_number; i++) {
			sorted_params[i] = unifrnd() / (double)error->model_parameter_number;
		}
		qsort(sorted_params, error->model_parameter_number, sizeof(double), comp_double);
		if (error->b_symmetric) {
			int j = 0;
			for (i = error->model_parameter_number - 1; i >= 0; i--) {
				parameters[j++] = sorted_params[i];
			}
		}
		else {
			int j = error->model_parameter_number - 1;
			parameters[error->model_parameter_diff] = sorted_params[j--];
			for (i = 1; i <= error->model_parameter_diff; i++) {
				parameters[error->model_parameter_diff - i] = sorted_params[j--];
				parameters[error->model_parameter_diff + i] = sorted_params[j--];
			}
		}

		pfm = fminsearch_new_with_eq(__loglikelihood_pairs_from_true_measure, error->model_parameter_number, error);
		pfm->tolx = 1e-9;
		pfm->tolf = 1e-9;
		fminsearch_min(pfm, parameters);
		double *re = fminsearch_get_minX(pfm);
		for (i = 0; i < error->model_parameter_number; i++) parameters[i] = re[i];
		log << std::endl << "Misclassification Matrix Search Result: (" << pfm->iters << " iterations)\n";
		log << "Score: " << *pfm->fv << std::endl;

		if (runs > 0) {
			if (!isnan(*pfm->fv) && !isinf(*pfm->fv) && abs(minscore - (*pfm->fv)) < pfm->tolf) {
				converged = 1;
			}
		}
		if (pfm->iters < pfm->maxiters) {
			if (*pfm->fv < minscore) {
				minscore = *pfm->fv;
				memcpy(bestrun_parameters, parameters, (error->model_parameter_number) * sizeof(double));
			}
			runs++;
		}
		fminsearch_free(pfm);
	} while (!converged && runs<max_runs);

	if (converged) {
		log << "score converged in " << runs << " runs." << std::endl;
	}
	else {
		log << "score failed to converge in " << max_runs << " runs." << std::endl;
		log << "best score: " << minscore << std::endl;
	}

	memory_free(parameters);
	error->estimates = bestrun_parameters;

	//memory_free(error);           // we are going to return these values
	memory_free_2dim((void**)observed_pairs, maxFamilySize + 1, maxFamilySize + 1, NULL);
	memory_free(sizeFreq);
	return error;
}



