#include "RcppExtended.hpp"
#include "quasarclTools.hpp"



//[[Rcpp::export]]
SEXP cppConvolve(SEXP quasarclPtr_, SEXP inputMatrix_, SEXP sizesVector_, SEXP filterVector_, SEXP outputMatrix_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix inputMatrix(inputMatrix_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	Rcpp::IntegerVector filterVector(filterVector_);
	Rcpp::NumericMatrix outputMatrix(outputMatrix_);
	
	const size_t width = inputMatrix.cols();
	const size_t height = inputMatrix.rows();
	const size_t N = width * height;
	const size_t M = sizesVector.size();
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferInput = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, M * sizeof(uint));
	cl::Buffer bufferFilter = cl::Buffer(context, CL_MEM_READ_ONLY, M * sizeof(double));
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	
	cl::copy(queue, inputMatrix.begin(), inputMatrix.end(), bufferInput);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
	cl::copy(queue, filterVector.begin(), filterVector.end(), bufferFilter);
	cl::copy(queue, outputMatrix.begin(), outputMatrix.end(), bufferOutput);
	
	quasarcl::convolve(*quasarclPtr, bufferInput, width, height, bufferSizes, bufferFilter, M, bufferOutput); 
	
	cl::copy(queue, bufferOutput, outputMatrix.begin(), outputMatrix.end());
	return outputMatrix;
}



//[[Rcpp::export]]
SEXP cppCopyIfNotInf(SEXP quasarclPtr_, SEXP inputMatrix_, SEXP filteredSize_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix inputMatrix(inputMatrix_);
	const size_t filteredSize = Rcpp::as<unsigned int>(filteredSize_);
	
	auto spectrumsSize = inputMatrix.cols();
	auto spectrumsNumber = inputMatrix.rows();
	const size_t N = spectrumsSize * spectrumsNumber;

	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferInput = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, filteredSize * spectrumsNumber * sizeof(double));
	
	cl::copy(queue, inputMatrix.begin(), inputMatrix.end(), bufferInput);
	
	quasarcl::copyIfNotInf(*quasarclPtr, bufferInput, spectrumsNumber, spectrumsSize, bufferOutput, filteredSize);
	
	Rcpp::NumericMatrix outputMatrix(spectrumsNumber, filteredSize);
	cl::copy(queue, bufferOutput, outputMatrix.begin(), outputMatrix.end());

	return outputMatrix;	
}



//[[Rcpp::export]]
SEXP cppCountIfNotInf(SEXP quasarclPtr_, SEXP inputMatrix_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix inputMatrix(inputMatrix_);

	auto spectrumsSize = inputMatrix.cols();
	auto spectrumsNumber = inputMatrix.rows();
	const size_t N = spectrumsSize * spectrumsNumber;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferInput = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, spectrumsNumber * sizeof(uint));

	cl::copy(queue, inputMatrix.begin(), inputMatrix.end(), bufferInput);
	
	quasarcl::countIfNotInf(*quasarclPtr, bufferInput, spectrumsNumber, spectrumsSize, bufferOutput); 
	
	Rcpp::IntegerVector outputVector(spectrumsNumber);
	cl::copy(queue, bufferOutput, outputVector.begin(), outputVector.end());
	return outputVector;
}



//[[Rcpp::export]]
SEXP cppReglin(SEXP quasarclPtr_, SEXP xMatrix_, SEXP yMatrix_, SEXP sizesVector_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix xMatrix(xMatrix_);
	Rcpp::NumericMatrix yMatrix(yMatrix_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	
	auto spectrumsNumber = xMatrix.rows();
	auto spectrumsSize = xMatrix.cols();
	const size_t N = spectrumsNumber * spectrumsSize;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferX = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferY = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(uint));
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, spectrumsNumber * sizeof(cl_double8));
	
	cl::copy(queue, xMatrix.begin(), xMatrix.end(), bufferX);
	cl::copy(queue, yMatrix.begin(), yMatrix.end(), bufferY);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
	
	
	quasarcl::reglin(*quasarclPtr, bufferX, bufferY, spectrumsNumber, spectrumsSize, bufferSizes, bufferOutput); 
	
	std::vector<cl_double8> outputVector(spectrumsNumber);
	cl::copy(queue, bufferOutput, outputVector.begin(), outputVector.end());
	return Rcpp::wrap(outputVector);
}



//[[Rcpp::export]]
SEXP cppReglinYax(SEXP quasarclPtr_, SEXP xMatrix_, SEXP yMatrix_, SEXP sizesVector_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix xMatrix(xMatrix_);
	Rcpp::NumericMatrix yMatrix(yMatrix_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	
	auto spectrumsNumber = xMatrix.rows();
	auto spectrumsSize = xMatrix.cols();
	const size_t N = spectrumsNumber * spectrumsSize;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferX = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferY = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(uint));
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, spectrumsNumber * sizeof(double));
	
	cl::copy(queue, xMatrix.begin(), xMatrix.end(), bufferX);
	cl::copy(queue, yMatrix.begin(), yMatrix.end(), bufferY);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
	
	quasarcl::reglin(*quasarclPtr, bufferX, bufferY, spectrumsNumber, spectrumsSize, bufferSizes, bufferOutput); 
	
	Rcpp::NumericVector outputVector(spectrumsNumber);
	cl::copy(queue, bufferOutput, outputVector.begin(), outputVector.end());
	return outputVector;
}



//[[Rcpp::export]]
SEXP cppChisq(SEXP quasarclPtr_, SEXP fMatrix_, SEXP yMatrix_, SEXP errorsMatrix_, SEXP sizesVector_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix fMatrix(fMatrix_);
	Rcpp::NumericMatrix yMatrix(yMatrix_);
	Rcpp::NumericMatrix errorsMatrix(errorsMatrix_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	
	auto spectrumsNumber = fMatrix.rows();
	auto spectrumsSize = fMatrix.cols();
	const size_t N = spectrumsNumber * spectrumsSize;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferF = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferY = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferErrors = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(uint));
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, spectrumsNumber * sizeof(double));
	
	cl::copy(queue, fMatrix.begin(), fMatrix.end(), bufferF);
	cl::copy(queue, yMatrix.begin(), yMatrix.end(), bufferY);
	cl::copy(queue, errorsMatrix.begin(), errorsMatrix.end(), bufferErrors);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
		
	quasarcl::chisq(*quasarclPtr, bufferF, bufferY, bufferErrors, spectrumsNumber, spectrumsSize, bufferSizes, bufferOutput); 
	
	Rcpp::NumericVector outputVector(spectrumsNumber);
	cl::copy(queue, bufferOutput, outputVector.begin(), outputVector.end());
	return outputVector;
}



//[[Rcpp::export]]
SEXP cppTrapz(SEXP quasarclPtr_, SEXP yMatrix_, SEXP xMatrix_, SEXP sizesVector_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix yMatrix(yMatrix_);
	Rcpp::NumericMatrix xMatrix(xMatrix_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	
	auto spectrumsNumber = xMatrix.rows();
	auto spectrumsSize = xMatrix.cols();
	const size_t N = spectrumsNumber * spectrumsSize;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferY = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferX = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(uint));
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, spectrumsNumber * sizeof(double));
	
	cl::copy(queue, yMatrix.begin(), yMatrix.end(), bufferY);
	cl::copy(queue, xMatrix.begin(), xMatrix.end(), bufferX);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
			
	quasarcl::trapz(*quasarclPtr, bufferY, bufferX, spectrumsNumber, spectrumsSize, bufferSizes, bufferOutput); 
	
	Rcpp::NumericVector outputVector(spectrumsNumber);
	cl::copy(queue, bufferOutput, outputVector.begin(), outputVector.end());
	return outputVector;
}

