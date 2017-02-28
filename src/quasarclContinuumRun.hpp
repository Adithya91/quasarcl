#ifndef QUASARCL_CONTINUUM_RUN_H_
#define QUASARCL_CONTINUUM_RUN_H_


#include "QuasarCL.hpp"

#include "quasarclBasics.hpp"
#include "quasarclTools.hpp"
#include "quasarclSpectrums.hpp"
#include "quasarclContinuum.hpp"


namespace quasarcl {
	

	struct ContinuumResults
	{
			cl::Buffer dcontinuumsMatrix;
			cl::Buffer continuumsMatrix;
			Rcpp::NumericVector chisqs;
			std::vector<cl_double8> c_reglin_result;
			std::vector<cl_double8> reglin_result;
	};

		
	
	inline ContinuumResults continuum(Rcpp::XPtr<QuasarCL> quasarclPtr, 
									  cl::Buffer& spectrumsMatrix, cl::Buffer& wavelengthsMatrix, 
									  cl::Buffer& errorsMatrix, cl::Buffer& sizes, const size_t size, 
									  const size_t spec_size, std::vector<cl_double2> windows,
									  cl_double ampWavelength)
	{
		auto queue = quasarclPtr->getQueue();
		auto context = quasarclPtr->getContext();
		
		//w original kopia danych z wavelengthsMatrix
		cl::Buffer wavelengthsMatrix_original = cl::Buffer(context, CL_MEM_READ_WRITE, ASTRO_OBJ_SPEC_SIZE * size * sizeof(cl_double));
		queue.enqueueCopyBuffer(wavelengthsMatrix, wavelengthsMatrix_original, 0, 0, ASTRO_OBJ_SPEC_SIZE * size * sizeof(cl_double));
		
		cl::Buffer windowsBuffer = cl::Buffer(context, CL_MEM_READ_ONLY, windows.size() * sizeof(cl_double2));
		cl::copy(queue, windows.begin(), windows.end(), windowsBuffer);		
		
		//Filtrowanie danych
		filterWithWavelengthWindows(*quasarclPtr, spectrumsMatrix, wavelengthsMatrix, errorsMatrix, sizes, size, windowsBuffer, windows.size());
		// ___________________________________
		
		filterNonpositive(*quasarclPtr, spectrumsMatrix, wavelengthsMatrix, errorsMatrix, sizes, size);
		
		cl::Buffer sizes_filtered = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(cl_uint));
		countIfNotInf(*quasarclPtr, spectrumsMatrix, size, ASTRO_OBJ_SPEC_SIZE, sizes_filtered);
	
		std::vector<cl_uint> sizes_filtered_vec(size, 0);
		cl::copy(queue, sizes_filtered, sizes_filtered_vec.begin(), sizes_filtered_vec.end());
		auto max = *std::max_element(sizes_filtered_vec.begin(), sizes_filtered_vec.end());

		size_t max_spec_filtered_size = max > 0 ? max : 64;
		max_spec_filtered_size = calcGlobalSize(64, max_spec_filtered_size);

		// Bufory na przefiltrowana dane.
		cl::Buffer spectrumsMatrix_filtered = cl::Buffer(context, CL_MEM_READ_WRITE, max_spec_filtered_size * size * sizeof(cl_double));
		cl::Buffer wavelengthsMatrix_filtered = cl::Buffer(context, CL_MEM_READ_WRITE, max_spec_filtered_size * size * sizeof(cl_double));

		copyIfNotInf(*quasarclPtr, spectrumsMatrix, size, ASTRO_OBJ_SPEC_SIZE, spectrumsMatrix_filtered, max_spec_filtered_size); //
		copyIfNotInf(*quasarclPtr, wavelengthsMatrix, size, ASTRO_OBJ_SPEC_SIZE, wavelengthsMatrix_filtered , max_spec_filtered_size);

		// Wolne bufory: spectrumsMatrix, wavelengthsMatrix

		auto errorsMatrix_filtered = spectrumsMatrix;
		copyIfNotInf(*quasarclPtr, errorsMatrix, size, ASTRO_OBJ_SPEC_SIZE, errorsMatrix_filtered, max_spec_filtered_size);

		//
		// Wolne bufory: wavelengthsMatrix, errorsMatrix

		auto spectrumsMatrix_filtered_copy = errorsMatrix;
		auto wavelengthsMatrix_filtered_copy = wavelengthsMatrix;

		// Zapisanie kopii buforów _filtered
		size_t byteSize = sizeof(cl_double) * max_spec_filtered_size * size;
		queue.enqueueCopyBuffer(spectrumsMatrix_filtered, spectrumsMatrix_filtered_copy, 0, 0, byteSize);
		queue.enqueueCopyBuffer(wavelengthsMatrix_filtered, wavelengthsMatrix_filtered_copy, 0, 0, byteSize);

		// ________________________________________________________________________ //
		// __wyznacznie "parametrów continuum" cpar________________________________ //
		// ________________________________________________________________________ //

		log10(*quasarclPtr, spectrumsMatrix_filtered_copy, size, max_spec_filtered_size);
		log10(*quasarclPtr, wavelengthsMatrix_filtered_copy, size, max_spec_filtered_size);
		
		// bufor na wyniki dopasowania prostej regresji (czyli jej współczynniki)
		// dla przefiltrowanych danych.
		cl::Buffer c_reglin_results = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(cl_double8));
		reglin(*quasarclPtr, wavelengthsMatrix_filtered_copy, spectrumsMatrix_filtered_copy, size, max_spec_filtered_size, sizes_filtered, c_reglin_results);

		
		if(ampWavelength > std::pow(DBL_MIN, 10.0001)) // Odrzucam takie lamp, że wynik log10(lamp) będzie równy -INFINITY
		{
			cl_double lampLog10 = static_cast<cl_double>(std::log10(static_cast<double>(ampWavelength)));
			minus(*quasarclPtr, wavelengthsMatrix_filtered_copy, size, max_spec_filtered_size, lampLog10);
		}
		
		// bufor na wyniki dopasowania prostej regresji (czyli jej współczynniki)
		// dla przefiltrowanych danych.
		cl::Buffer reglin_results = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(cl_double8));
		reglin(*quasarclPtr, wavelengthsMatrix_filtered_copy, spectrumsMatrix_filtered_copy, size, max_spec_filtered_size, sizes_filtered, reglin_results);
		
		// Wolne bufory: errorsMatrix, wavelengthsMatrix.

		// Przeliczenie parametrów i błędów w przestrzeni widma.
		// Modyfikacja współczynnika b w c_reglin_results: b = 10^b.
		fixReglinResults(*quasarclPtr, c_reglin_results, reglin_results, size);
		
		//
		auto continuumsMatrix_filtered = wavelengthsMatrix_filtered; // Zapisuje w tym samym buforze, z którego biorę dane (kernel na to pozwala).		
		calcCw(*quasarclPtr, wavelengthsMatrix_filtered, continuumsMatrix_filtered, max_spec_filtered_size, size, c_reglin_results);
		
		auto reducedChisqs_filtered = errorsMatrix; // To co jest w errors nie będzie więcej potrzebne.
		chisq(*quasarclPtr, spectrumsMatrix_filtered, continuumsMatrix_filtered, errorsMatrix_filtered, size, max_spec_filtered_size, sizes_filtered, reducedChisqs_filtered);
		reduceContinuumChisqs(*quasarclPtr, reducedChisqs_filtered, sizes_filtered, size);

		// ________________________________________________________________________ //
		// ___obliczanie continuum i błędów continuum______________________________ //
		// ________________________________________________________________________ //

		//
		// errorsMatrix_filtered niepotrzebne.
		// Wolne: wavelengthsMatrix, spectrumsMatrix

		//
		auto dcontinuumsMatrix = wavelengthsMatrix; // To co jest w lambdas nie będzie więcej potrzebne.
		auto continuumsMatrix = spectrumsMatrix; // Jak wyżej.

		calcCfunDcfun(*quasarclPtr, wavelengthsMatrix_original, dcontinuumsMatrix, continuumsMatrix, size, c_reglin_results, reglin_results);

		// ________________________________________________________________________ //
		// ___kopiowanie wyników na host___________________________________________ //
		// ________________________________________________________________________ //

		Rcpp::NumericVector reducedChisqs_filteredVec(size);
		std::vector<cl_double8> reglin_resultsVec(size);
		std::vector<cl_double8> c_reglin_resultsVec(size);

		cl::copy(queue, reducedChisqs_filtered, reducedChisqs_filteredVec.begin(), reducedChisqs_filteredVec.end());
		cl::copy(queue, reglin_results, reglin_resultsVec.begin(), reglin_resultsVec.end());
		cl::copy(queue, c_reglin_results, c_reglin_resultsVec.begin(), c_reglin_resultsVec.end());
		
		ContinuumResults result = {
			dcontinuumsMatrix,
			continuumsMatrix,
			reducedChisqs_filteredVec,
			reglin_resultsVec,
			c_reglin_resultsVec
		};
		return result;
	}

}

#endif