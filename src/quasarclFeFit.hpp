#ifndef QUASARCL_FE_FIT_RUN_H_
#define QUASARCL_FE_FIT_RUN_H_


#include "QuasarCL.hpp"

#include "quasarclTools.hpp"
#include "quasarclSpectrums.hpp"
#include "quasarclFeWindows.hpp"
#include "quasarclBasics.hpp"

using namespace Rcpp;

namespace quasarcl {
	

	struct FitResults
	{
			cl::Buffer feTemplateMatrix;
			Rcpp::NumericVector scaleRates;
			Rcpp::IntegerVector sizes_fewindows;
			Rcpp::NumericVector reducedChisqs_fewindows;
			Rcpp::NumericVector reducedChisqs_full;
			Rcpp::NumericVector ews_full;
			Rcpp::NumericVector reducedChisqs_feRange;
			Rcpp::NumericVector ews_feRange;
	};
	
		
	
	struct Data
	{		
			const Rcpp::NumericVector& spectrumsMatrixHost;
			const Rcpp::NumericVector& wavelengthsMatrixHost;
			const Rcpp::NumericVector& errorsMatrixHost;
			const Rcpp::NumericVector& continuumsMatrixHost;
			const Rcpp::IntegerVector& sizes;
	};

	

	struct Buffers
	{
			cl::Buffer& spectrumsMatrix;
			cl::Buffer& wavelengthsMatrix;
			cl::Buffer& errorsMatrix;
			cl::Buffer& continuumsMatrix;
			cl::Buffer& sizes;
	};		


	
	inline cl::Buffer calcFeTemplateScaleRates(Rcpp::XPtr<QuasarCL> quasarclPtr, cl::Buffer spectrumsMatrix, cl::Buffer templateFeMatrix,
			cl::Buffer sizes, size_t size, size_t max_spectrum_size,
			Rcpp::List fitParameters)
	{
		// obliczanie współczynników a dla prostej regresji liniowej y = a * x
		// prosta dla dopasowania szablonu żelaza do widma
		auto queue = quasarclPtr->getQueue();
		auto context = quasarclPtr->getContext();
		
		cl::Buffer reglinYaxResults = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(cl_double));
		reglinYax(*quasarclPtr, templateFeMatrix, spectrumsMatrix, size, max_spectrum_size,
				sizes, reglinYaxResults);

		// Wymnożenie współczynników z reglinYaxResults przez stała wartość skalującą feScaleRate
		// podaną w parametrach algorytmu
		Rcpp::NumericVector reglinYaxResultsVec(size);
		cl::copy(queue, reglinYaxResults, reglinYaxResultsVec.begin(), reglinYaxResultsVec.end());
		cl_double feScaleRate = Rcpp::as<cl_double>(fitParameters["feScaleRate"]);

		std::transform(reglinYaxResultsVec.begin(), reglinYaxResultsVec.end(), reglinYaxResultsVec.begin(),
				[&](cl_double x)
				{	return x * feScaleRate;});

		cl::copy(queue, reglinYaxResultsVec.begin(), reglinYaxResultsVec.end(), reglinYaxResults);

		return reglinYaxResults;
	}




	inline cl::Buffer calcReducedChisqs(Rcpp::XPtr<QuasarCL> quasarclPtr, cl::Buffer& fs, cl::Buffer& ys, cl::Buffer errors,
			cl::Buffer sizes, size_t size, size_t max_spectrum_size)
	{
		auto queue = quasarclPtr->getQueue();
		auto context = quasarclPtr->getContext();
		
		// Obliczenie chi kwadrat
		cl::Buffer chisqResults = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(cl_double));
		chisq(*quasarclPtr, fs, ys, errors, size, max_spectrum_size, sizes, chisqResults);
		// Podzielenie chi^2 przez ilość elementów
		// w celu uzyskanie ZREDUKOWANEGO CHI KWADRAT
		reduceFeChisqs(*quasarclPtr, chisqResults, sizes, size);
		return chisqResults;
	}




	inline Rcpp::NumericVector cpuConvolve(Rcpp::NumericVector signal, Rcpp::NumericVector kernel, bool same)
	{
		Rcpp::NumericVector result(signal.size() + kernel.size() - 1, 0);

		for (size_t n = 0; n < result.size(); n++)
		{
			size_t kmin, kmax, k;

			kmin = (n >= kernel.size() - 1) ? n - (kernel.size() - 1) : 0;
			kmax = (n < signal.size() - 1) ? n : signal.size() - 1;

			for (k = kmin; k <= kmax; k++)
			{
				result[n] += signal[k] * kernel[n - k];
			}
		}

		if (same)
		{
			size_t kernel_center = kernel.size() % 2 == 0 ? (kernel.size() - 1) / 2 : kernel.size() / 2;
			result.erase(result.begin(), result.begin() + kernel_center);
			result.erase(result.end() - (kernel.size() - 1 - kernel_center), result.end());
		}

		return result;
	}
	
	
	
	inline cl::Buffer calcFeTemplateMatrix(Rcpp::XPtr<QuasarCL> quasarclPtr, cl::Buffer wavelengthsMatrix, cl::Buffer sizes, size_t size,
			Rcpp::List feTemplate, Rcpp::List fitParameters)
	{
		auto queue = quasarclPtr->getQueue();
		auto context = quasarclPtr->getContext();
		
		Rcpp::NumericVector templateFeValuesVec;
		if (Rcpp::as<double>(fitParameters["fwhmn"]) > Rcpp::as<double>(fitParameters["fwhmt"]))
		{
			// Prędkość światła [m/s]
			const double C = 299792458.0;

			double sincov = std::pow(std::pow( Rcpp::as<double>(fitParameters["fwhmn"]), 2.0) - std::pow( Rcpp::as<double>(fitParameters["fwhmt"]), 2.0), 0.5)
					/ 2.0;
			sincov /= pow((2.0 * std::log(2.0)), 0.5) * C / 1e3;

			Rcpp::NumericVector feGauss = Rcpp::clone(Rcpp::as<Rcpp::NumericVector>(feTemplate["wavelengths"]));

			std::transform(feGauss.begin(), feGauss.end(), feGauss.begin(),
					[](cl_double x) -> cl_double
					{
						return std::log10(x);
					});
				
			// pik gaussa z offsetem
			cl_double a = 1.0 / sincov / std::pow(2 * M_PI, 0.5);
			cl_double b = feGauss[feGauss.size() / 2];
			cl_double c = sincov;
			cl_double d = 0.0;
			std::transform(feGauss.begin(), feGauss.end(), feGauss.begin(),
					[&](cl_double x) -> cl_double
					{
						return a * std::exp(-0.5 * (x-b) * (x-b)/pow(c, 2.0)) + d;
					});
			
			Rcpp::NumericVector feValues = feTemplate["values"];
			templateFeValuesVec = cpuConvolve(feValues, feGauss, true);
		}
		else
		{
			templateFeValuesVec = feTemplate["values"];
		}
		cl::Buffer templateFeMatrix = cl::Buffer(context, CL_MEM_READ_WRITE, size * ASTRO_OBJ_SPEC_SIZE * sizeof(cl_double));
		{
			cl::Buffer zeros = cl::Buffer(context, CL_MEM_READ_WRITE, size * ASTRO_OBJ_SPEC_SIZE * sizeof(cl_double));
			{
				Rcpp::NumericVector zerosVec(size * ASTRO_OBJ_SPEC_SIZE, 0.0);
				cl::copy(queue, zerosVec.begin(), zerosVec.end(), zeros);
			}
			
			Rcpp::NumericVector templateFeWavelengthsVec = feTemplate["wavelengths"];
			
			cl::Buffer templateFeValues = cl::Buffer(context, CL_MEM_READ_ONLY, templateFeValuesVec.size() * sizeof(cl_double));
			cl::Buffer templateFeLambda = cl::Buffer(context, CL_MEM_READ_ONLY, templateFeWavelengthsVec.size() * sizeof(cl_double));

			cl::copy(queue, templateFeValuesVec.begin(), templateFeValuesVec.end(), templateFeValues);		
			cl::copy(queue, templateFeWavelengthsVec.begin(), templateFeWavelengthsVec.end(), templateFeLambda);

			// Dodawanie szablonu żelaza do zer z interpolacją zgodnie z długościami fal widm.
			addSpectrum(*quasarclPtr, zeros, wavelengthsMatrix, sizes, size,
					templateFeLambda, templateFeValues, templateFeValuesVec.size(),
					templateFeMatrix);
		}

		return templateFeMatrix;
	}
	
	
	
	inline FitResults feFit(Rcpp::XPtr<QuasarCL> quasarclPtr, Data spectrumsData, Buffers spectrumsBuffers, 
							const size_t size, Rcpp::List feTemplate, std::vector<cl_double2> feWindows,
							Rcpp::List fitParameters)
	{
		auto queue = quasarclPtr->getQueue();
		auto context = quasarclPtr->getContext();
		
		// ________________________________________________________________________ //
		// __bufory ze spectrumsBuffer_____________________________________________ //
		// ________________________________________________________________________ //

		cl::Buffer spectrumsMatrix = spectrumsBuffers.spectrumsMatrix;
		cl::Buffer wavelengthsMatrix = spectrumsBuffers.wavelengthsMatrix;
		cl::Buffer errorsMatrix = spectrumsBuffers.errorsMatrix;
		cl::Buffer continuumsMatrix = spectrumsBuffers.continuumsMatrix;
		cl::Buffer sizes = spectrumsBuffers.sizes;

		// ________________________________________________________________________ //
		// __przygotowanie/przetwarzanie szablonu żelaza___________________________ //
		// ________________________________________________________________________ //

		// Macierz, w której w każdej kolumnie będzie znajdował się szablon żelaza
		// "dopasowany" do długości fal dla każdego wima.
		//
		// Można powiedzieć, że wycinamy ten kawałek szablonu, który jest w zakresie długości wal widma
		// i jeszcze odpowiednio go dopasowujemy do długości fal widma (interpolacja).
		//
		auto templateFeMatrix = calcFeTemplateMatrix(quasarclPtr, wavelengthsMatrix, sizes, size, feTemplate, fitParameters);
		// ________________________________________________________________________ //
		// __zachowanie kopii obliczonego szablonu żelaza__________________________ //
		// ________________________________________________________________________ //
		cl::Buffer templateFeMatrixCopy = cl::Buffer(context, CL_MEM_READ_WRITE, size * ASTRO_OBJ_SPEC_SIZE * sizeof(cl_double));
		queue.enqueueCopyBuffer(templateFeMatrix, templateFeMatrixCopy, 0, 0, size * ASTRO_OBJ_SPEC_SIZE * sizeof(cl_double));

		// ________________________________________________________________________ //
		// __jeżeli continuum nie zostało odjęte od widm to odejmujemy_____________ //
		// ________________________________________________________________________ //	
		
		if (!Rcpp::as<bool>(fitParameters["isSubC"]))
		{
			minus(*quasarclPtr, spectrumsMatrix, size, ASTRO_OBJ_SPEC_SIZE, continuumsMatrix, spectrumsMatrix);
		}
		
		// ________________________________________________________________________ //
		// __oblicznia odpowiednie dla wybranego typu obszaru______________________ //
		// ________________________________________________________________________ //

		std::string type = Rcpp::as<std::string>(fitParameters["fitType"]);
		
		if (type == "WIN" || type == "FWIN")
		{
			// cl::Buffer na okna
			cl::Buffer feWindowsBuffer = cl::Buffer(context, CL_MEM_READ_ONLY, feWindows.size() * sizeof(cl_double2));
			cl::copy(queue, feWindows.begin(), feWindows.end(), feWindowsBuffer);

			// Filtrowanie danych po oknach (długości fal)

			// Widmo, długości fal, błędy.
			//
			
			filterWithWavelengthWindows(*quasarclPtr, spectrumsMatrix, wavelengthsMatrix, errorsMatrix, sizes, size,
										feWindowsBuffer, feWindows.size());

			
			// Żelazo i kontinuum
			//
			// specs jest już przefiltrowane, więc w "złych" miejscach jest INIFNITY
			filterInfs(*quasarclPtr, spectrumsMatrix, continuumsMatrix, templateFeMatrix, sizes, size);
		}
		if (type == "FWIN")
		{
			// Filtrowanie zer, usuwanie danych z tych długości fal
			// dla których wartość szablonu żelaza (tj. z odpowiadającej kolumny templateFeMatrix) jest równa zero.

			// Widmo, długości fal, błędy.
			//
			filterZeros(*quasarclPtr, templateFeMatrix, wavelengthsMatrix, errorsMatrix, sizes, size);

			// Żelazo i kontinuum
			//
			// specs jest już przefiltrowane, więc w "złych" miejscach jest INIFNITY
			filterInfs(*quasarclPtr, templateFeMatrix, continuumsMatrix, spectrumsMatrix, sizes, size);
		}

		// ________________________________________________________________________ //
		// __Skracanie danych po filtracji_________________________________________ //
		// ________________________________________________________________________ //

		// Jeżeli filtrowaliśmy dane, no to warto działać tylko na tych które przeszły
		// filtrowanie.
		auto spectrumsMatrix_filtered = spectrumsMatrix;
		auto wavelengthsMatrix_filtered = wavelengthsMatrix;
		auto errorsMatrix_filtered = errorsMatrix;
		auto continuumsMatrix_filtered = continuumsMatrix;
		auto templateFeMatrix_filtered = templateFeMatrix;
		auto sizes_fewindows_filtered = sizes;
		size_t max_spectrum_filtered_size = ASTRO_OBJ_SPEC_SIZE;

		if (type == "WIN" || type == "FWIN")
		{
			// Rozmiary widm po filtrowaniu.
			sizes_fewindows_filtered = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(cl_uint));
			// Zliczanie elementów po przefiltrowaniu.
			countIfNotInf(*quasarclPtr, spectrumsMatrix, size, ASTRO_OBJ_SPEC_SIZE, sizes_fewindows_filtered);
			
			
			// Pobranie rozmiarów i wybranie największego (max), następnie
			// obliczenie najbliżej wielokrotności 64 większej od maksimum.
			// Potrzebne, żeby adresy elementów początkowych były wielokrotnością 64;
			Rcpp::IntegerVector sizes_filtered_vec(size, 0);
			cl::copy(queue, sizes_fewindows_filtered, sizes_filtered_vec.begin(), sizes_filtered_vec.end());
			auto max = *std::max_element(sizes_filtered_vec.begin(), sizes_filtered_vec.end());
			
			size_t max_spectrum_filtered_size = max > 0 ? max : 64;			
			max_spectrum_filtered_size = calcGlobalSize(64, max_spectrum_filtered_size);			

			copyIfNotInf(*quasarclPtr, spectrumsMatrix, size, ASTRO_OBJ_SPEC_SIZE, spectrumsMatrix, max_spectrum_filtered_size);
			copyIfNotInf(*quasarclPtr, wavelengthsMatrix, size, ASTRO_OBJ_SPEC_SIZE, wavelengthsMatrix, max_spectrum_filtered_size);
			copyIfNotInf(*quasarclPtr, errorsMatrix, size, ASTRO_OBJ_SPEC_SIZE, errorsMatrix, max_spectrum_filtered_size);
			copyIfNotInf(*quasarclPtr, continuumsMatrix, size, ASTRO_OBJ_SPEC_SIZE, continuumsMatrix, max_spectrum_filtered_size);
			copyIfNotInf(*quasarclPtr, templateFeMatrix, size, ASTRO_OBJ_SPEC_SIZE, templateFeMatrix, max_spectrum_filtered_size);
			

			spectrumsMatrix_filtered = spectrumsMatrix;
			wavelengthsMatrix_filtered = wavelengthsMatrix;
			errorsMatrix_filtered = errorsMatrix;
			continuumsMatrix_filtered = continuumsMatrix;
			templateFeMatrix_filtered = templateFeMatrix;
		}

		// ________________________________________________________________________ //
		// __obliczanie współczynników skalujących dla szablonów żelaza____________ //
		// ________________________________________________________________________ //
		//std::cout<<"Before calculating scale rates"<<std::endl;
		cl::Buffer scaleRates = calcFeTemplateScaleRates(quasarclPtr, spectrumsMatrix_filtered, templateFeMatrix_filtered,
															   sizes_fewindows_filtered, size, max_spectrum_filtered_size, fitParameters);


		// _________________________________________________________________________ //
		// __wymnożenie szablonu żelaza przez współczynniki skalujące_________________ //
		// ___________________________________________________________________________ //

		// przefiltrowanego
		multiplyCol(*quasarclPtr, templateFeMatrix_filtered, size, max_spectrum_filtered_size, scaleRates,
					templateFeMatrix_filtered);

		// pełnego szablonu żelaza
		multiplyCol(*quasarclPtr, templateFeMatrixCopy, size, ASTRO_OBJ_SPEC_SIZE, scaleRates, templateFeMatrixCopy);

		// ________________________________________________________________________ //
		// __jakość dopasowania dla danych z okien (przefiltrowanych)______________ //
		// ________________________________________________________________________ //

		cl::Buffer reducedChisqs_filtered = calcReducedChisqs(quasarclPtr, spectrumsMatrix_filtered, templateFeMatrix_filtered, errorsMatrix_filtered,
															  sizes_fewindows_filtered, size, max_spectrum_filtered_size);

		// ________________________________________________________________________ //
		// __załadowanie oryginalnych danych potrzebnych do dalszych obliczeń______ //
		// ________________________________________________________________________ //

		// Załadowanie pełnych wersji widm, długości fal, błędów i kontinuum
		cl::copy(queue, spectrumsData.spectrumsMatrixHost.begin(), spectrumsData.spectrumsMatrixHost.end(),
				spectrumsMatrix);
		cl::copy(queue, spectrumsData.wavelengthsMatrixHost.begin(), spectrumsData.wavelengthsMatrixHost.end(),
				wavelengthsMatrix);
		cl::copy(queue, spectrumsData.errorsMatrixHost.begin(), spectrumsData.errorsMatrixHost.end(),
				errorsMatrix);
		cl::copy(queue, spectrumsData.continuumsMatrixHost.begin(), spectrumsData.continuumsMatrixHost.end(),
				continuumsMatrix);

		// jeżeli continuum nie zostało odjęte od widm to odejmujemy
		if (!Rcpp::as<bool>(fitParameters["isSubC"]))
		{
			minus(*quasarclPtr, spectrumsMatrix, size, ASTRO_OBJ_SPEC_SIZE, continuumsMatrix, spectrumsMatrix);
		}

		// ________________________________________________________________________ //
		// __jakość dopasowania dla pełnego widma__________________________________ //
		// ________________________________________________________________________ //

		cl::Buffer reducedChisqs_full = calcReducedChisqs(quasarclPtr, spectrumsMatrix, templateFeMatrixCopy, errorsMatrix,
														  sizes, size, ASTRO_OBJ_SPEC_SIZE);

		// ________________________________________________________________________ //
		// __pole powierzchni dla pełnego szablonu żelaza__________________________ //
		// ________________________________________________________________________ //
		cl::Buffer ews_full = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(cl_double));
		{
			// templateFeMatrix już nie jest potrzebne, więc można go wykorzystać.
			auto temp = templateFeMatrix;
			divide(*quasarclPtr, templateFeMatrixCopy, size, ASTRO_OBJ_SPEC_SIZE, continuumsMatrix, temp);
			trapz(*quasarclPtr, temp, wavelengthsMatrix, size, ASTRO_OBJ_SPEC_SIZE, sizes, ews_full);
		}

		// ________________________________________________________________________ //
		// __przywrócenie pełnej wersji templateFeMatrix z kopii___________________ //
		// ________________________________________________________________________ //

		queue.enqueueCopyBuffer(templateFeMatrixCopy, templateFeMatrix, 0, 0, size * ASTRO_OBJ_SPEC_SIZE * sizeof(cl_double));

		// ________________________________________________________________________ //
		// __jakość dopasowania i ew (pole pow) dla głównego pasma żelaza__________ //
		// ________________________________________________________________________ //

		cl::Buffer reducedChisqs_feRange;
		cl::Buffer ews_feRange;


		if (type == "WIN" || type == "FWIN")
		{
			// ________________________________________________________________________ //
			// __fltrowanie danych, tylko główne pasmo żelaza__________________________ //
			// ________________________________________________________________________ //
			cl::Buffer sizes_feRange_filtered = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(cl_uint));

			cl_double2 feFitRange = Rcpp::as<cl_double2>(fitParameters["feFitRange"]);
			std::vector<cl_double2> feFitRangeVec(1, feFitRange);

			cl::Buffer feFitRangeBuffer = cl::Buffer(context, CL_MEM_READ_ONLY, feFitRangeVec.size() * sizeof(cl_double2));
			cl::copy(queue, feFitRangeVec.begin(), feFitRangeVec.end(), feFitRangeBuffer);

			filterWithWavelengthWindows(*quasarclPtr, spectrumsMatrix, wavelengthsMatrix, errorsMatrix, sizes, size,
										feFitRangeBuffer, feFitRangeVec.size());
			filterInfs(*quasarclPtr, spectrumsMatrix, continuumsMatrix, templateFeMatrix, sizes, size);

			// Zliczanie elementów po przefiltrowaniu.
			countIfNotInf(*quasarclPtr, spectrumsMatrix, size, ASTRO_OBJ_SPEC_SIZE, sizes_feRange_filtered);

			// ________________________________________________________________________ //
			// __Skracanie danych po filtracji_________________________________________ //
			// ________________________________________________________________________ //

			// Pobranie rozmiarów i wybranie największego (max), następnie
			// obliczenie najbliżej wielokrotności 64 większej od maksimum.
			// Potrzebne, żeby adresy elementów początkowych były wielokrotnością 64;
			Rcpp::IntegerVector sizes_filtered_vec(size, 0);
			cl::copy(queue, sizes_feRange_filtered, sizes_filtered_vec.begin(), sizes_filtered_vec.end());
			auto max = *std::max_element(sizes_filtered_vec.begin(), sizes_filtered_vec.end());

			max_spectrum_filtered_size = max > 0 ? max : 64;
			max_spectrum_filtered_size = calcGlobalSize(64, max_spectrum_filtered_size);

			copyIfNotInf(*quasarclPtr, spectrumsMatrix, size, ASTRO_OBJ_SPEC_SIZE, spectrumsMatrix, max_spectrum_filtered_size);
			copyIfNotInf(*quasarclPtr, wavelengthsMatrix, size, ASTRO_OBJ_SPEC_SIZE, wavelengthsMatrix, max_spectrum_filtered_size);
			copyIfNotInf(*quasarclPtr, errorsMatrix, size, ASTRO_OBJ_SPEC_SIZE, errorsMatrix, max_spectrum_filtered_size);
			copyIfNotInf(*quasarclPtr, continuumsMatrix, size, ASTRO_OBJ_SPEC_SIZE, continuumsMatrix, max_spectrum_filtered_size);
			copyIfNotInf(*quasarclPtr, templateFeMatrix, size, ASTRO_OBJ_SPEC_SIZE, templateFeMatrix, max_spectrum_filtered_size);

			// ________________________________________________________________________ //
			// __jakość dopasowania dla danych z głównego pasma żelaza_________________ //
			// ________________________________________________________________________ //

			reducedChisqs_feRange = calcReducedChisqs(quasarclPtr, spectrumsMatrix, templateFeMatrix, errorsMatrix,
													  sizes_feRange_filtered, size, max_spectrum_filtered_size);

			// ________________________________________________________________________ //
			// __pole powierzchni dla pełnego szablonu żelaza__________________________ //
			// ________________________________________________________________________ //

			{
				ews_feRange = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(cl_double));
				divide(*quasarclPtr, templateFeMatrix, size, ASTRO_OBJ_SPEC_SIZE, continuumsMatrix, templateFeMatrix);
				trapz(*quasarclPtr, templateFeMatrix, wavelengthsMatrix, size, ASTRO_OBJ_SPEC_SIZE, sizes_feRange_filtered, ews_feRange);
			}
		}
		// fitParameters.fitType == Type::FULL
		else
		{
			reducedChisqs_feRange = reducedChisqs_full;
			ews_feRange = ews_full;
		}
		
		// ________________________________________________________________________ //
		// __kopiowanie wyników do wektorów________________________________________ //
		// ________________________________________________________________________ //


		Rcpp::NumericVector scaleRatesVec(size);
		Rcpp::IntegerVector sizes_fewindows_filteredVec(size);
		Rcpp::NumericVector reducedChisqs_filteredVec(size);
		Rcpp::NumericVector reducedChisqs_fullVec(size);
		Rcpp::NumericVector ews_fullVec(size);
		Rcpp::NumericVector reducedChisqs_feRangeVec(size);
		Rcpp::NumericVector ews_feRangeVec(size);

		cl::copy(queue, scaleRates, scaleRatesVec.begin(), scaleRatesVec.end());
		cl::copy(queue, sizes_fewindows_filtered, sizes_fewindows_filteredVec.begin(), sizes_fewindows_filteredVec.end());
		cl::copy(queue, reducedChisqs_filtered, reducedChisqs_filteredVec.begin(), reducedChisqs_filteredVec.end());
		cl::copy(queue, reducedChisqs_full, reducedChisqs_fullVec.begin(), reducedChisqs_fullVec.end());
		cl::copy(queue, ews_full, ews_fullVec.begin(), ews_fullVec.end());
		cl::copy(queue, reducedChisqs_feRange, reducedChisqs_feRangeVec.begin(), reducedChisqs_feRangeVec.end());
		cl::copy(queue, ews_feRange, ews_feRangeVec.begin(), ews_feRangeVec.end());
		
		FitResults results = {
			templateFeMatrixCopy,
			scaleRatesVec,
			sizes_fewindows_filteredVec,
			reducedChisqs_filteredVec,
			reducedChisqs_fullVec,
			ews_fullVec,
			reducedChisqs_feRangeVec,
			ews_feRangeVec
		};
		return results;
	}

}

#endif