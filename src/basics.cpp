#include "quasarclBasics.hpp"


//[[Rcpp::export]]
SEXP cppLog10(SEXP quasarclPtr_, SEXP inputMatrix_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix inputMatrix(inputMatrix_);
	
	const size_t width = inputMatrix.cols();   //długość widma
	const size_t height = inputMatrix.rows();  //liczba kwazarów
	const size_t N = width * height;

	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferInput = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::copy(queue, inputMatrix.begin(), inputMatrix.end(), bufferInput);
	
	quasarcl::log10(*quasarclPtr, bufferInput, height, width);
	
	Rcpp::NumericMatrix outputMatrix(height, width);
	cl::copy(queue, bufferInput, outputMatrix.begin(), outputMatrix.end());
	return outputMatrix;
}


//[[Rcpp::export]]
SEXP cppMinusMatrix(SEXP quasarclPtr_, SEXP inputMatrix_, SEXP subtrahendMatrix_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix inputMatrix(inputMatrix_);
	Rcpp::NumericMatrix subtrahendMatrix(subtrahendMatrix_);

	const size_t width = inputMatrix.cols();   //długość widma
	const size_t height = inputMatrix.rows();  //liczba kwazarów
	const size_t N = width * height;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferInput = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferSubtrahend = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	
	cl::copy(queue, inputMatrix.begin(), inputMatrix.end(), bufferInput);
	cl::copy(queue, subtrahendMatrix.begin(), subtrahendMatrix.end(), bufferSubtrahend);
	
	quasarcl::minus(*quasarclPtr, bufferInput, height, width, bufferSubtrahend, bufferOutput);
	
	Rcpp::NumericMatrix outputMatrix(height, width);
	cl::copy(queue, bufferOutput, outputMatrix.begin(), outputMatrix.end());
	return outputMatrix;
}


//[[Rcpp::export]]
SEXP cppMinusScalar(SEXP quasarclPtr_, SEXP inputMatrix_, SEXP scalar_) {

	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix inputMatrix(inputMatrix_);
	double scalar = Rcpp::as<double>(scalar_);
	
	const size_t width = inputMatrix.cols();   //długość widma
	const size_t height = inputMatrix.rows();  //liczba kwazarów
	const size_t N = width * height;
	
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferInput = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::copy(queue, inputMatrix.begin(), inputMatrix.end(), bufferInput);
	
	quasarcl::minus(*quasarclPtr, bufferInput, height, width, scalar);
		
	Rcpp::NumericMatrix outputMatrix = Rcpp::NumericMatrix(height, width);
	cl::copy(queue, bufferInput, outputMatrix.begin(), outputMatrix.end());
	return outputMatrix;
}


//[[Rcpp::export]]
SEXP cppMultiplyCol(SEXP quasarclPtr_, SEXP inputMatrix_, SEXP vector_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix inputMatrix(inputMatrix_);
	Rcpp::NumericVector vector(vector_);
	
	const size_t width = inputMatrix.cols();   //długość widma
	const size_t height = inputMatrix.rows();  //liczba kwazarów
	const size_t N = width * height;
	const size_t M = vector.size();

	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferInput = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferVector = cl::Buffer(context, CL_MEM_READ_ONLY, M * sizeof(double));
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	
	cl::copy(queue, inputMatrix.begin(), inputMatrix.end(), bufferInput);
	cl::copy(queue, vector.begin(), vector.end(), bufferVector);
	
	quasarcl::multiplyCol(*quasarclPtr, bufferInput, height, width, bufferVector, bufferOutput);
	
	Rcpp::NumericMatrix outputMatrix(height, width);
	cl::copy(queue, bufferOutput, outputMatrix.begin(), outputMatrix.end());
	return outputMatrix;
}



//[[Rcpp::export]]
SEXP cppTranspose(SEXP quasarclPtr_, SEXP inputMatrix_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix inputMatrix(inputMatrix_);
	
	const size_t width = inputMatrix.cols();  //liczba kwazarów
	const size_t height = inputMatrix.rows(); //długość widma
	const size_t N = width * height;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferInput = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));

	cl::copy(queue, inputMatrix.begin(), inputMatrix.end(), bufferInput);
	
	quasarcl::transpose(*quasarclPtr, bufferInput, height, width, bufferOutput);
	
	Rcpp::NumericMatrix outputMatrix(width, height);
	cl::copy(queue, bufferOutput, outputMatrix.begin(), outputMatrix.end());	
	return outputMatrix;
}



//[[Rcpp::export]]
SEXP cppDivide(SEXP quasarclPtr_, SEXP inputMatrix_, SEXP divisorMatrix_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix inputMatrix(inputMatrix_);
	Rcpp::NumericMatrix divisorMatrix(divisorMatrix_);
	
	const size_t width = inputMatrix.cols();   //długość widma
	const size_t height = inputMatrix.rows();  //liczba kwazarów
	const size_t N = width * height;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferInput = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferDivisor = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	
	cl::copy(queue, inputMatrix.begin(), inputMatrix.end(), bufferInput);
	cl::copy(queue, divisorMatrix.begin(), divisorMatrix.end(), bufferDivisor);
	
	quasarcl::divide(*quasarclPtr, bufferInput, height, width, bufferDivisor, bufferOutput);
	
	Rcpp::NumericMatrix outputMatrix(height, width);
	cl::copy(queue, bufferOutput, outputMatrix.begin(), outputMatrix.end());
	return outputMatrix;
}
	
