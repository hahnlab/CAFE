#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <stdexcept>
#include <iterator>
#include "simerror.h"

int backup_original_count(pCafeFamily pcf)
{
	int idx = 0;

	if (pcf->countbackup == NULL) {
		pcf->countbackup = (int **)memory_new(pcf->flist->size, sizeof(int*));
	}
	for (idx = 0; idx < pcf->flist->size; idx++)
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[idx];
		pcf->countbackup[idx] = pitem->count;
	}
	return 0;
}


int restore_original_count(pCafeFamily pcf)
{
	int idx = 0;

	if (((pCafeFamilyItem)pcf->flist->array[0])->count != pcf->countbackup[0]) {
		for (idx = 0; idx < pcf->flist->size; idx++)
		{
			pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[idx];
			memory_free(pitem->count);
			pitem->count = pcf->countbackup[idx];
		}
		memory_free(pcf->countbackup);
		pcf->countbackup = NULL;
	}
	return 0;
}

/// random sampling based on misclassification probability.
/// This may have a bug in it as there is nothing inherently preventing r from going over the limit
/// Is there an unstated assumption that that misclassification_probabilities add up to 1?
size_t get_random(std::vector<double> misclassification_probability)
{
	size_t r;
	double cumul = 0;
	double rnd = unifrnd();
	for (r = 0; r < misclassification_probability.size(); r++)
	{
		cumul += misclassification_probability[r];
		if (rnd <= cumul) break;
	}

	return r;

}

void write_strings(std::ostream& ost, char** items, int size, std::string delimiter)
{
	std::copy(items, items+size - 1, std::ostream_iterator<std::string>(ost, delimiter.c_str()));
	ost << items[size-1];
}

void write_species_counts(pCafeFamily pcf, std::ostream& ost)
{
	ost << "Desc\tFamily ID\t";
	write_strings(ost, pcf->species, pcf->num_species, "\t");
	ost << "\n";

	for (int idx = 0; idx < pcf->flist->size; idx++)
	{
		pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[idx];
		ost << pitem->desc << "\t" << pitem->id << "\t" << pitem->count[0];
		for (int n = 1; n < pcf->num_species; n++)
		{
			ost << "\t" << pitem->count[n];
		}
		ost << "\n";
	}
}

void simulate_misclassification(pCafeFamily pcf)
{
	for (int idx = 0; idx < pcf->flist->size; idx++)
	{
		int n;
		pCafeFamilyItem pitem = (pCafeFamilyItem)pcf->flist->array[idx];
		std::vector<int> originalcount(pcf->num_species);
		std::copy(pitem->count, pitem->count + pcf->num_species, originalcount.begin());

		pitem->count = (int*)memory_new(pcf->num_species, sizeof(int));
		for (n = 0; n < pcf->num_species; n++)
		{
			int j = originalcount[n];     // column idx
			pErrorStruct errormodel = pcf->error_ptr[n];        // right error model
			std::vector<double> misclassification_probability(errormodel->maxfamilysize);

			// copy the conditional misclassification probability given true size (column j of errormatrix)
			for (int i = 0; i < errormodel->maxfamilysize; i++) {
				misclassification_probability[i] = errormodel->errormatrix[i][j];
			}

			pitem->count[n] = get_random(misclassification_probability);
		}
	}
}

/**
to check how simex is working
generate data with additional error (Ystar) by applying the misclassification matrix to the true data
( random sampling of true factors=f, with probabilities following the column=f of the misclassification matrix  )
then run simex on the variable with error (Ystar), and compare it with the true model (based on Y) and the naive model(based on Ystar).


first find the naive estimator by assuming there is no error in the data.
then we will add even more error and store the estimates for each dataset with increasing error (k = 0.5, 1, 1.5, 2)
so in the end we have (1+number of lambda) * (number of parameters) estimates.
for each i = k we run the estimate j = B number of times.
each time add error by applying the misclassification matrix = errormatrix^k
update the estimates based on error added data,
store the mean estimates for all B runs.
then predict the estimates at k=-1

**/
double simerror(pCafeFamily pcf, std::string prefix, int repeat)
{
	for (int run = 0; run<repeat; run++)
	{
		backup_original_count(pcf);

		std::ostringstream trainfile;
		trainfile << prefix << "_" << run << ".erred";
		std::ofstream ofst(trainfile.str().c_str());
		if (!ofst)
		{
			throw std::runtime_error("ERROR(simerror): Cannot open " + trainfile.str() + " in write mode.");
		}

		simulate_misclassification(pcf);

		write_species_counts(pcf, ofst);

		restore_original_count(pcf);

	}
	return 0;
}

