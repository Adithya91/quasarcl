#include "RcppExtended.hpp"
#include "quasarclSpectrums.hpp"


//[[Rcpp::export]]
SEXP cppGenerateWavelengthsMatrix(SEXP quasarclPtr_, SEXP abzVector_, SEXP spectrumsSize_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	std::vector<cl_double4> abzVector = Rcpp::as<std::vector<cl_double4> >(abzVector_);
	const size_t spectrumsSize = Rcpp::as<unsigned int>(spectrumsSize_);

	auto spectrumsNumber = abzVector.size();
	const size_t N = spectrumsNumber * spectrumsSize;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferAbz = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(cl_double4));
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));

	cl::copy(queue, abzVector.begin(), abzVector.end(), bufferAbz);
	
	quasarcl::generateWavelengthsMatrix(*quasarclPtr, bufferAbz, spectrumsNumber, bufferOutput); 
	
	Rcpp::NumericMatrix outputMatrix(spectrumsSize, spectrumsNumber);
	cl::copy(queue, bufferOutput, outputMatrix.begin(), outputMatrix.end());
	return outputMatrix;
}



//[[Rcpp::export]]
SEXP cppAddSpectrum(SEXP quasarclPtr_, SEXP wavelengthsMatrix_, SEXP spectrumsMatrix_, SEXP sizesVector_, SEXP toAddWavelengthsVector_, SEXP toAddSpectrumVector_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix wavelengthsMatrix(wavelengthsMatrix_);
	Rcpp::NumericMatrix spectrumsMatrix(spectrumsMatrix_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	Rcpp::NumericVector toAddWavelengthsVector(toAddWavelengthsVector_);
	Rcpp::NumericVector toAddSpectrumVector(toAddSpectrumVector_);
	
	auto spectrumsNumber = wavelengthsMatrix.rows();
	auto spectrumsSize = wavelengthsMatrix.cols();
	const size_t N = spectrumsNumber * spectrumsSize;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferWavelengths = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferSpectrums = cl::Buffer(context, CL_MEM_READ_ONLY, N * sizeof(double));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(uint));
	cl::Buffer bufferToAddWavelenghts = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(double));
	cl::Buffer bufferToAddSpectrum = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(double));
	cl::Buffer bufferOutput = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	
	cl::copy(queue, wavelengthsMatrix.begin(), wavelengthsMatrix.end(), bufferWavelengths);
	cl::copy(queue, spectrumsMatrix.begin(), spectrumsMatrix.end(), bufferSpectrums);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
	cl::copy(queue, toAddWavelengthsVector.begin(), toAddWavelengthsVector.end(), bufferToAddWavelenghts);
	cl::copy(queue, toAddSpectrumVector.begin(), toAddSpectrumVector.end(), bufferToAddSpectrum);
	
	quasarcl::addSpectrum(*quasarclPtr, bufferWavelengths, bufferSpectrums, bufferSizes, spectrumsNumber, bufferToAddWavelenghts, bufferToAddSpectrum, spectrumsSize, bufferOutput); 
	
	Rcpp::NumericMatrix outputMatrix(spectrumsNumber, spectrumsSize);
	cl::copy(queue, bufferOutput, outputMatrix.begin(), outputMatrix.end());
	return outputMatrix;
}



//[[Rcpp::export]]
SEXP cppFilterWithWavelengthWindows(SEXP quasarclPtr_, SEXP spectrumsMatrix_, SEXP wavelengthsMatrix_, SEXP errorsMatrix_, SEXP sizesVector_, SEXP windowsVector_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix spectrumsMatrix(spectrumsMatrix_);
	Rcpp::NumericMatrix wavelengthsMatrix(wavelengthsMatrix_);
	Rcpp::NumericMatrix errorsMatrix(errorsMatrix_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	std::vector<cl_double2> windowsVector = Rcpp::as<std::vector<cl_double2> >(windowsVector_);
	
	auto spectrumsNumber = spectrumsMatrix.rows();
	auto spectrumsSize = spectrumsMatrix.cols();
	auto windowsNumber = windowsVector.size();
	const size_t N = spectrumsNumber * spectrumsSize;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferSpectrums = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferWavelengths = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferErrors = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(uint));
	cl::Buffer bufferWindows = cl::Buffer(context, CL_MEM_READ_ONLY, windowsNumber * sizeof(cl_double2));
	
	cl::copy(queue, wavelengthsMatrix.begin(), wavelengthsMatrix.end(), bufferWavelengths);
	cl::copy(queue, spectrumsMatrix.begin(), spectrumsMatrix.end(), bufferSpectrums);
	cl::copy(queue, errorsMatrix.begin(), errorsMatrix.end(), bufferErrors);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
	cl::copy(queue, windowsVector.begin(), windowsVector.end(), bufferWindows);
	
	quasarcl::filterWithWavelengthWindows(*quasarclPtr, bufferSpectrums, bufferWavelengths, bufferErrors, bufferSizes, spectrumsNumber, bufferWindows, windowsNumber); 
	
	Rcpp::NumericMatrix outputSpectrumsMatrix(spectrumsNumber, spectrumsSize);
	Rcpp::NumericMatrix outputWavelengthsMatrix(spectrumsNumber, spectrumsSize);
	Rcpp::NumericMatrix outputErrorsMatrix(spectrumsNumber, spectrumsSize);
	
	cl::copy(queue, bufferSpectrums, outputSpectrumsMatrix.begin(), outputSpectrumsMatrix.end());
	cl::copy(queue, bufferWavelengths, outputWavelengthsMatrix.begin(), outputWavelengthsMatrix.end());
	cl::copy(queue, bufferErrors, outputErrorsMatrix.begin(), outputErrorsMatrix.end());
	
	return Rcpp::List::create(Rcpp::Named("spectrumsMatrix") = outputSpectrumsMatrix,
							  Rcpp::Named("wavelengthsMatrix") = outputWavelengthsMatrix,
							  Rcpp::Named("errorsMatrix") = outputErrorsMatrix);
}



//[[Rcpp::export]]
SEXP cppFilterNonpositive(SEXP quasarclPtr_, SEXP spectrumsMatrix_, SEXP aMatrix_, SEXP bMatrix_, SEXP sizesVector_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix spectrumsMatrix(spectrumsMatrix_);
	Rcpp::NumericMatrix aMatrix(aMatrix_);
	Rcpp::NumericMatrix bMatrix(bMatrix_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	
	auto spectrumsNumber = spectrumsMatrix.rows();
	auto spectrumsSize = spectrumsMatrix.cols();
	const size_t N = spectrumsNumber * spectrumsSize;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferSpectrums = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferA = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferB = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(size_t));
	
	cl::copy(queue, spectrumsMatrix.begin(), spectrumsMatrix.end(), bufferSpectrums);
	cl::copy(queue, aMatrix.begin(), aMatrix.end(), bufferA);
	cl::copy(queue, bMatrix.begin(), bMatrix.end(), bufferB);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
	
	quasarcl::filterNonpositive(*quasarclPtr, bufferSpectrums, bufferA, bufferB, bufferSizes, spectrumsNumber); 
	
	Rcpp::NumericMatrix outputSpectrumsMatrix(spectrumsNumber, spectrumsSize);
	Rcpp::NumericMatrix outputAMatrix(spectrumsNumber, spectrumsSize);
	Rcpp::NumericMatrix outputBMatrix(spectrumsNumber, spectrumsSize);
	
	cl::copy(queue, bufferSpectrums, outputSpectrumsMatrix.begin(), outputSpectrumsMatrix.end());
	cl::copy(queue, bufferA, outputAMatrix.begin(), outputAMatrix.end());
	cl::copy(queue, bufferB, outputBMatrix.begin(), outputBMatrix.end());
	
	return Rcpp::List::create(Rcpp::Named("spectrumsMatrix") = outputSpectrumsMatrix,
							  Rcpp::Named("aMatrix") = outputAMatrix,
							  Rcpp::Named("bMatrix") = outputBMatrix);
}



//[[Rcpp::export]]
SEXP cppFilterZeros(SEXP quasarclPtr_, SEXP spectrumsMatrix_, SEXP aMatrix_, SEXP bMatrix_, SEXP sizesVector_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix spectrumsMatrix(spectrumsMatrix_);
	Rcpp::NumericMatrix aMatrix(aMatrix_);
	Rcpp::NumericMatrix bMatrix(bMatrix_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	
	auto spectrumsNumber = spectrumsMatrix.rows();
	auto spectrumsSize = spectrumsMatrix.cols();
	const size_t N = spectrumsNumber * spectrumsSize;

	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferSpectrums = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferA = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferB = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(uint));
	
	cl::copy(queue, spectrumsMatrix.begin(), spectrumsMatrix.end(), bufferSpectrums);
	cl::copy(queue, aMatrix.begin(), aMatrix.end(), bufferA);
	cl::copy(queue, bMatrix.begin(), bMatrix.end(), bufferB);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
	
	quasarcl::filterZeros(*quasarclPtr, bufferSpectrums, bufferA, bufferB, bufferSizes, spectrumsNumber);
	
	Rcpp::NumericMatrix outputSpectrumsMatrix(spectrumsNumber, spectrumsSize);
	Rcpp::NumericMatrix outputAMatrix(spectrumsNumber, spectrumsSize);
	Rcpp::NumericMatrix outputBMatrix(spectrumsNumber, spectrumsSize);
	
	cl::copy(queue, bufferSpectrums, outputSpectrumsMatrix.begin(), outputSpectrumsMatrix.end());
	cl::copy(queue, bufferA, outputAMatrix.begin(), outputAMatrix.end());
	cl::copy(queue, bufferB, outputBMatrix.begin(), outputBMatrix.end());
	
	return Rcpp::List::create(Rcpp::Named("spectrumsMatrix") = outputSpectrumsMatrix,
							  Rcpp::Named("aMatrix") = outputAMatrix,
							  Rcpp::Named("bMatrix") = outputBMatrix);
}



//[[Rcpp::export]]
SEXP cppFilterInfs(SEXP quasarclPtr_, SEXP spectrumsMatrix_, SEXP aMatrix_, SEXP bMatrix_, SEXP sizesVector_)
{
	Rcpp::XPtr<quasarcl::QuasarCL> quasarclPtr(quasarclPtr_);
	Rcpp::NumericMatrix spectrumsMatrix(spectrumsMatrix_);
	Rcpp::NumericMatrix aMatrix(aMatrix_);
	Rcpp::NumericMatrix bMatrix(bMatrix_);
	Rcpp::IntegerVector sizesVector(sizesVector_);
	
	auto spectrumsNumber = spectrumsMatrix.rows();
	auto spectrumsSize = spectrumsMatrix.cols();
	const size_t N = spectrumsNumber * spectrumsSize;
	
	auto context = quasarclPtr->getContext();
	auto queue = quasarclPtr->getQueue();
	
	cl::Buffer bufferSpectrums = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferA = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferB = cl::Buffer(context, CL_MEM_READ_WRITE, N * sizeof(double));
	cl::Buffer bufferSizes = cl::Buffer(context, CL_MEM_READ_ONLY, spectrumsNumber * sizeof(uint));
	
	cl::copy(queue, spectrumsMatrix.begin(), spectrumsMatrix.end(), bufferSpectrums);
	cl::copy(queue, aMatrix.begin(), aMatrix.end(), bufferA);
	cl::copy(queue, bMatrix.begin(), bMatrix.end(), bufferB);
	cl::copy(queue, sizesVector.begin(), sizesVector.end(), bufferSizes);
			
	quasarcl::filterInfs(*quasarclPtr, bufferSpectrums, bufferA, bufferB, bufferSizes, spectrumsNumber); 
	
	Rcpp::NumericMatrix outputSpectrumsMatrix(spectrumsNumber, spectrumsSize);
	Rcpp::NumericMatrix outputAMatrix(spectrumsNumber, spectrumsSize);
	Rcpp::NumericMatrix outputBMatrix(spectrumsNumber, spectrumsSize);
	
	cl::copy(queue, bufferSpectrums, outputSpectrumsMatrix.begin(), outputSpectrumsMatrix.end());
	cl::copy(queue, bufferA, outputAMatrix.begin(), outputAMatrix.end());
	cl::copy(queue, bufferB, outputBMatrix.begin(), outputBMatrix.end());
	
	return Rcpp::List::create(Rcpp::Named("spectrumsMatrix") = outputSpectrumsMatrix,
							  Rcpp::Named("aMatrix") = outputAMatrix,
							  Rcpp::Named("bMatrix") = outputBMatrix);
}
