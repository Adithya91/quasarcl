#include "RcppExtended.hpp"
#include "quasarclParameterization.hpp"
	
	
	
//[[Rcpp::export]]	
SEXP cppParameterization(SEXP quasarclPtr_, SEXP spectrumsMatrix_, SEXP errorsMatrix_, 
						 SEXP sizesVector_, SEXP abzVector_, SEXP spectralLinesList_, 
						 SEXP continuumWindowsVector_, SEXP ampWavelength_, 
						 SEXP feWindowsVector_, SEXP feTemplateList_, SEXP fitParametersList_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix spectrumsMatrix(spectrumsMatrix_);
	Rcpp::NumericMatrix errorsMatrix(errorsMatrix_);
	Rcpp::IntegerVector sizesVector(sizesVector_);

	std::vector<cl_double4> abzVector = Rcpp::as<std::vector<cl_double4>>(abzVector_);
	Rcpp::List feTemplate(feTemplateList_);
	
	std::vector<cl_double2> continuumWindowsVector = Rcpp::as<std::vector<cl_double2>>(continuumWindowsVector_);
	std::vector<cl_double2> feWindowsVector = Rcpp::as<std::vector<cl_double2>>(feWindowsVector_);
	Rcpp::List spectralLines(spectralLinesList_);
	double ampWavelength = Rcpp::as<double>(ampWavelength_);
	Rcpp::List fitParameters(fitParametersList_);

	
	auto results = quasarcl::parameterization(quasarclPtr, spectrumsMatrix, errorsMatrix, sizesVector,
											 abzVector, spectralLines, continuumWindowsVector,
											 ampWavelength, feWindowsVector, feTemplate, 
											 fitParameters); 
	return results;
}



//[[Rcpp::export]]
SEXP cppFitElement(SEXP quasarclPtr_, SEXP specLinesMatrix_, 
				   SEXP continuumsMatrix_, SEXP wavelengthsMatrix_,
				   SEXP errorsMatrix_, SEXP sizesVector_, SEXP element_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix specLinesMatrix(specLinesMatrix_);
	Rcpp::NumericMatrix continuumsMatrix(continuumsMatrix_);
	Rcpp::NumericMatrix wavelengthsMatrix(wavelengthsMatrix_);
	Rcpp::NumericMatrix errorsMatrix(errorsMatrix_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	
	Rcpp::List element(element_);
		
	auto spectrumsNumber = wavelengthsMatrix.rows();
	auto spectrumsSize = wavelengthsMatrix.cols();
	const size_t N = spectrumsNumber * spectrumsSize;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferSpecLines = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferSpecLinesCopy = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferContinuums = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferWavelengths = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferErrors = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_WRITE, spectrumsNumber * sizeof(uint));
	
	cl::copy(queue, specLinesMatrix.begin(), specLinesMatrix.end(), bufferSpecLines);
	cl::copy(queue, continuumsMatrix.begin(), continuumsMatrix.end(), bufferContinuums);
	cl::copy(queue, wavelengthsMatrix.begin(), wavelengthsMatrix.end(), bufferWavelengths);
	cl::copy(queue, errorsMatrix.begin(), errorsMatrix.end(), bufferErrors);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
	
	queue.enqueueCopyBuffer(bufferSpecLinesCopy, bufferSpecLines, 0, 0, N * sizeof(double));
	
	return quasarcl::fitElement(quasarclPtr, bufferSpecLines, bufferSpecLinesCopy,
								bufferContinuums, bufferWavelengths, wavelengthsMatrix,
								bufferErrors, bufferSizes, sizesVector, spectrumsNumber, element); 
}