#include "RcppExtended.hpp"
#include "quasarclContinuum.hpp"
	
	
	
	//[[Rcpp::export]]
	SEXP cppFixReglinResults(SEXP quasarclPtr_, SEXP cReglinResultsVector_, SEXP reglinResultsVector_)
	{
		Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
		std::vector<cl_double8> cReglinResultsVector = Rcpp::as<std::vector<cl_double8> >(cReglinResultsVector_);
		std::vector<cl_double8> reglinResultsVector = Rcpp::as<std::vector<cl_double8> >(reglinResultsVector_);
		const size_t size = cReglinResultsVector.size();
		
		auto context = quasarclPtr->getContext();
		auto queue = quasarclPtr->getQueue();

		cl::Buffer bufferCReglinResults = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(cl_double8));
		cl::Buffer bufferReglinResults = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(cl_double8));

		cl::copy(queue, cReglinResultsVector.begin(), cReglinResultsVector.end(), bufferCReglinResults);
		cl::copy(queue, reglinResultsVector.begin(), reglinResultsVector.end(), bufferReglinResults);		
		
		quasarcl::fixReglinResults(*quasarclPtr, bufferCReglinResults, bufferReglinResults, size);
		
		std::vector<cl_double8> outputCReglinResultsVector(size);
		std::vector<cl_double8> outputReglinResultsVector(size);
		
		cl::copy(queue, bufferCReglinResults, outputCReglinResultsVector.begin(), outputCReglinResultsVector.end());
		cl::copy(queue, bufferReglinResults, outputReglinResultsVector.begin(), outputReglinResultsVector.end());
			
		return Rcpp::List::create(Rcpp::Named("cReglinResults") = Rcpp::wrap(outputCReglinResultsVector),
								  Rcpp::Named("reglinResults") = Rcpp::wrap(outputReglinResultsVector));
		
	}
	
	

	//[[Rcpp::export]]
	SEXP cppCalcCfunDcfun(SEXP quasarclPtr_, SEXP wavelengthsMatrix_, SEXP cReglinResultsVector_, SEXP reglinResultsVector_)
	{
		Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
		Rcpp::NumericMatrix wavelengthsMatrix(wavelengthsMatrix_);
		
		std::vector<cl_double8> cReglinResultsVector = Rcpp::as<std::vector<cl_double8>>(cReglinResultsVector_);
		std::vector<cl_double8> reglinResultsVector = Rcpp::as<std::vector<cl_double8>>(reglinResultsVector_);
		
		const size_t spectrumsSize = wavelengthsMatrix.cols();
		const size_t spectrumsNumber = wavelengthsMatrix.rows();
		const size_t N = spectrumsSize * spectrumsNumber;
		
		auto context = quasarclPtr->getContext();
		auto queue = quasarclPtr->getQueue();
		
		cl::Buffer bufferWavelengths = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(cl_double));
		cl::Buffer bufferDcontinuums = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(cl_double));
		cl::Buffer bufferContinuums = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(cl_double));
		cl::Buffer bufferCReglinResults = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(cl_double8));
		cl::Buffer bufferReglinResults = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(cl_double8));
		
		cl::copy(queue, wavelengthsMatrix.begin(), wavelengthsMatrix.end(), bufferWavelengths);
		cl::copy(queue, cReglinResultsVector.begin(), cReglinResultsVector.end(), bufferCReglinResults);
		cl::copy(queue, reglinResultsVector.begin(), reglinResultsVector.end(), bufferReglinResults);
				
		quasarcl::calcCfunDcfun(*quasarclPtr, bufferWavelengths, bufferDcontinuums, bufferContinuums, spectrumsNumber, 
								bufferCReglinResults, bufferReglinResults);
		
		Rcpp::NumericMatrix outputDContinuumsMatrix(spectrumsNumber, spectrumsSize);
		Rcpp::NumericMatrix outputContinuumsMatrix(spectrumsNumber, spectrumsSize);
		
		cl::copy(queue, bufferDcontinuums, outputDContinuumsMatrix.begin(), outputDContinuumsMatrix.end());
		cl::copy(queue, bufferContinuums, outputContinuumsMatrix.begin(), outputContinuumsMatrix.end());
				
		return Rcpp::List::create(Rcpp::Named("dcontinuums") = outputDContinuumsMatrix, 
								  Rcpp::Named("continuums") = outputContinuumsMatrix);
	}
	
	
	
	//[[Rcpp::export]]
	SEXP cppCalcCw(SEXP quasarclPtr_, SEXP wavelengthsMatrix_, SEXP cReglinResultsVector_)
	{
		Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
		Rcpp::NumericMatrix wavelengthsMatrix(wavelengthsMatrix_);
		
		std::vector<cl_double8> cReglinResultsVector = Rcpp::as<std::vector<cl_double8> >(cReglinResultsVector_);
		
		auto spectrumsNumber = wavelengthsMatrix.rows();
		auto spectrumsSize = wavelengthsMatrix.cols();
		const size_t N = spectrumsNumber * spectrumsSize;
		
		auto context = quasarclPtr->getContext();
		auto queue = quasarclPtr->getQueue();
		
		cl::Buffer bufferWavelengths = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(cl_double));
		cl::Buffer bufferContinuums = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(cl_double));
		cl::Buffer bufferCReglinResults = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(cl_double8));
		
		cl::copy(queue, wavelengthsMatrix.begin(), wavelengthsMatrix.end(), bufferWavelengths);
		cl::copy(queue, cReglinResultsVector.begin(), cReglinResultsVector.end(), bufferCReglinResults);
			
		quasarcl::calcCw(*quasarclPtr, bufferWavelengths, bufferContinuums, spectrumsSize, spectrumsNumber, 
						 bufferCReglinResults);
		
		Rcpp::NumericMatrix outputContinuumsMatrix(spectrumsNumber, spectrumsSize);
		cl::copy(queue, bufferContinuums, outputContinuumsMatrix.begin(), outputContinuumsMatrix.end());		
		return outputContinuumsMatrix;
	} 
	
	
	
	//[[Rcpp::export]]
	SEXP cppReduceContinuumChisqs(SEXP quasarclPtr_, SEXP chisqsVector_, SEXP sizesVector_)
	{
		Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
		Rcpp::NumericVector chisqsVector(chisqsVector_);
		Rcpp::IntegerVector sizesVector(sizesVector_);
		
		const size_t size = chisqsVector.size();
		
		auto context = quasarclPtr->getContext();
		auto queue = quasarclPtr->getQueue();
		
		cl::Buffer bufferChisqs = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(cl_double));
		cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, size * sizeof(cl_uint));
		
		cl::copy(queue, chisqsVector.begin(), chisqsVector.end(), bufferChisqs);
		cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
		
		quasarcl::reduceContinuumChisqs(*quasarclPtr, bufferChisqs, bufferSizes, size);
		
		Rcpp::NumericVector outputChisqsVector(size);
		cl::copy(queue, bufferChisqs, outputChisqsVector.begin(), outputChisqsVector.end());
		return outputChisqsVector;
	}  