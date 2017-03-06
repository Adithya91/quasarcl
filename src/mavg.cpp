#include "quasarclMAVG.hpp"


//[[Rcpp::export]]
SEXP cppSimpleMAVG(SEXP quasarclPtr_, SEXP inputMatrix_, SEXP windowWidth_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix inputMatrix(inputMatrix_);
	const unsigned int windowWidth = Rcpp::as<unsigned int>(windowWidth_);
	
	auto spectrumsNumber = inputMatrix.rows();
	auto spectrumsSize = inputMatrix.cols();
	const size_t N = spectrumsNumber * spectrumsSize;

	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferInput = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(cl_double));;
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(cl_double));
	
	cl::copy(queue, inputMatrix.begin(), inputMatrix.end(), bufferInput);
	
	quasarcl::simpleMAVG(*quasarclPtr, bufferInput, spectrumsNumber, spectrumsSize, bufferOutput, windowWidth); 
	
	Rcpp::NumericMatrix outputMatrix(spectrumsNumber, spectrumsSize);
	cl::copy(queue, bufferOutput, outputMatrix.begin(), outputMatrix.end());
	return outputMatrix;
}



//[[Rcpp::export]]
SEXP cppCenteredMAVG(SEXP quasarclPtr_, SEXP inputMatrix_, SEXP colsSizesVector_, SEXP windowWidth_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix inputMatrix(inputMatrix_);
	Rcpp::IntegerVector colsSizesVector(colsSizesVector_);
	
	
	const unsigned int windowWidth = Rcpp::as<unsigned int>(windowWidth_);
	
	auto spectrumsNumber = inputMatrix.rows();
	auto spectrumsSize = inputMatrix.cols();
	const size_t N = spectrumsNumber * spectrumsSize;
	

	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferInput = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(cl_double));
	cl::Buffer bufferColsSizes = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(cl_uint));
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(cl_double));
	
	cl::copy(queue, inputMatrix.begin(), inputMatrix.end(), bufferInput);
	cl::copy(queue, colsSizesVector.begin(), colsSizesVector.end(), bufferColsSizes);
	
	quasarcl::centeredMAVG(*quasarclPtr, bufferInput, spectrumsNumber, spectrumsSize, bufferColsSizes, bufferOutput, windowWidth);  
	
	Rcpp::NumericMatrix outputMatrix(spectrumsNumber, spectrumsSize);
	cl::copy(queue, bufferOutput, outputMatrix.begin(), outputMatrix.end());
	return outputMatrix;
}

	
