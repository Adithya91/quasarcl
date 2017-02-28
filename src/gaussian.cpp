#include "RcppExtended.hpp"
#include "quasarclGaussian.hpp"

	
//[[Rcpp::export]]
SEXP cppFitGaussian(SEXP quasarclPtr_, SEXP yMatrix_, SEXP xMatrix_, SEXP sizesVector_)
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
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, spectrumsNumber * sizeof(cl_double4));
	
	cl::copy(queue, xMatrix.begin(), xMatrix.end(), bufferX);
	cl::copy(queue, yMatrix.begin(), yMatrix.end(), bufferY);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
	
	quasarcl::fitGaussian(*quasarclPtr, bufferY, bufferX, bufferSizes, spectrumsNumber, bufferOutput);
	
	std::vector<cl_double4> outputVector(spectrumsNumber);	
	cl::copy(queue, bufferOutput, outputVector.begin(), outputVector.end());
	return Rcpp::wrap(outputVector);
	
}



//[[Rcpp::export]]
SEXP cppCalcGaussian(SEXP quasarclPtr_, SEXP xMatrix_, SEXP gaussianParamsVector_, SEXP sizesVector_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix xMatrix(xMatrix_);
	std::vector<cl_double4> gaussianParamsVector = Rcpp::as<std::vector<cl_double4> >(gaussianParamsVector_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	
	auto spectrumsNumber = xMatrix.rows();
	auto spectrumsSize = xMatrix.cols();
	
	const size_t N = spectrumsNumber * spectrumsSize;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferX = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferGaussianParams = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(cl_double4));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(uint));
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	
	cl::copy(queue, xMatrix.begin(), xMatrix.end(), bufferX);
	cl::copy(queue, gaussianParamsVector.begin(), gaussianParamsVector.end(), bufferGaussianParams);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
	
	quasarcl::calcGaussian(*quasarclPtr, bufferX, bufferGaussianParams,  bufferSizes, spectrumsNumber, spectrumsSize, bufferOutput);
	
	Rcpp::NumericMatrix outputMatrix(spectrumsNumber, spectrumsSize);	
	cl::copy(queue, bufferOutput, outputMatrix.begin(), outputMatrix.end());
	return outputMatrix;
	
}



//[[Rcpp::export]]
SEXP cppCalcGaussianChisqs(SEXP quasarclPtr_, SEXP xMatrix_, SEXP yMatrix_, SEXP errorsMatrix_, SEXP gaussianParamsVector_, SEXP sizesVector_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix xMatrix(xMatrix_);
	Rcpp::NumericMatrix yMatrix(yMatrix_);
	Rcpp::NumericMatrix errorsMatrix(errorsMatrix_);
	std::vector<cl_double4> gaussianParamsVector = Rcpp::as<std::vector<cl_double4> >(gaussianParamsVector_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	
	auto spectrumsNumber = xMatrix.rows();
	auto spectrumsSize = xMatrix.cols();
	
	const size_t N = spectrumsNumber * spectrumsSize;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferX = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferY = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferErrors = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferGaussianParams = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(cl_double4));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(uint));
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, spectrumsNumber * sizeof(double));
	
	cl::copy(queue, xMatrix.begin(), xMatrix.end(), bufferX);
	cl::copy(queue, yMatrix.begin(), yMatrix.end(), bufferY);
	cl::copy(queue, errorsMatrix.begin(), errorsMatrix.end(), bufferErrors);
	cl::copy(queue, gaussianParamsVector.begin(), gaussianParamsVector.end(), bufferGaussianParams);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
	
	quasarcl::calcGaussianChisq(*quasarclPtr, bufferX, bufferY, bufferErrors, bufferGaussianParams, bufferSizes, spectrumsNumber, bufferOutput);
	
	Rcpp::NumericVector outputVector(spectrumsNumber);	
	cl::copy(queue, bufferOutput, outputVector.begin(), outputVector.end());
	return outputVector;
	
}



//[[Rcpp::export]]
SEXP cppCalcGaussianFWHM(SEXP quasarclPtr_, SEXP gaussianParamsVector_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	std::vector<cl_double4> gaussianParamsVector = Rcpp::as<std::vector<cl_double4> >(gaussianParamsVector_);
	
	const size_t size = gaussianParamsVector.size();
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferGaussianParams = cl::Buffer(context, CL_MEM_READ_ONLY, size * sizeof(cl_double4));
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, size * sizeof(double));
	
	cl::copy(queue, gaussianParamsVector.begin(), gaussianParamsVector.end(), bufferGaussianParams);
	
	quasarcl::calcGaussianFWHM(*quasarclPtr, bufferGaussianParams, bufferOutput, size);
	
	Rcpp::NumericVector outputVector(size);	
	cl::copy(queue, bufferOutput, outputVector.begin(), outputVector.end());
	return outputVector;
	
}