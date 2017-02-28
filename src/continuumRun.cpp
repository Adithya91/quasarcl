#include "RcppExtended.hpp"
#include "quasarclContinuumRun.hpp"



	//[[Rcpp::export]]
	SEXP cppContinuum(SEXP quasarclPtr_, SEXP spectrumsMatrix_, SEXP wavelengthsMatrix_, SEXP errorsMatrix_, 
					  SEXP sizesVector_, SEXP windowsVector_, SEXP ampWavelength_)
	{
		Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
		Rcpp::NumericMatrix spectrumsMatrix(spectrumsMatrix_);
		Rcpp::NumericMatrix wavelengthsMatrix(wavelengthsMatrix_);
		Rcpp::NumericMatrix errorsMatrix(errorsMatrix_);
		Rcpp::IntegerVector sizesVector(sizesVector_);
		
		std::vector<cl_double2> windowsVector = Rcpp::as<std::vector<cl_double2> >(windowsVector_);
		double ampWavelength = Rcpp::as<double>(ampWavelength_);
		
		auto spectrumsNumber = spectrumsMatrix.rows();
		auto spectrumsSize = spectrumsMatrix.cols();
		const size_t N = spectrumsNumber * spectrumsSize;
		
		auto context = quasarclPtr->getContext();
		auto queue = quasarclPtr->getQueue();
		
		cl::Buffer bufferSpectrums = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
		cl::Buffer bufferWavelengths = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
		cl::Buffer bufferErrors = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
		cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_WRITE, spectrumsNumber * sizeof(uint));
		
		cl::copy(queue, spectrumsMatrix.begin(), spectrumsMatrix.end(), bufferSpectrums);
		cl::copy(queue, wavelengthsMatrix.begin(), wavelengthsMatrix.end(), bufferWavelengths);
		cl::copy(queue, errorsMatrix.begin(), errorsMatrix.end(), bufferErrors);
		cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
		
		auto results = quasarcl::continuum(quasarclPtr, bufferSpectrums, bufferWavelengths, bufferErrors, bufferSizes, 
										   spectrumsNumber, spectrumsSize, windowsVector, ampWavelength); 
		
		Rcpp::NumericMatrix dcontinuumsMatrix(spectrumsNumber, spectrumsSize);
		Rcpp::NumericMatrix continuumsMatrix(spectrumsNumber, spectrumsSize);
		
		cl::copy(queue, results.dcontinuumsMatrix, dcontinuumsMatrix.begin(), dcontinuumsMatrix.end());
		cl::copy(queue, results.continuumsMatrix, continuumsMatrix.begin(), continuumsMatrix.end());
		
		return Rcpp::List::create(Rcpp::Named("dcontinuumsMatrix") = dcontinuumsMatrix,
								  Rcpp::Named("continuumsMatrix") = continuumsMatrix,
								  Rcpp::Named("chisqs") = results.chisqs,
								  Rcpp::Named("c_reglin_result") = Rcpp::wrap(results.c_reglin_result),
								  Rcpp::Named("reglin_result") = Rcpp::wrap(results.reglin_result));
	}