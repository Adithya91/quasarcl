#' @useDynLib quasar
#' @importFrom Rcpp evalCpp



#' @export
rFeFit <- function(quasarcl, spectrumsMatrix, wavelengthsMatrix, errorsMatrix, continuumsMatrix,
					  sizesVector, feTemplate, windows = NULL, fitParameters = NULL) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(dim(spectrumsMatrix), dim(wavelengthsMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(dim(spectrumsMatrix), dim(errorsMatrix)))
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(dim(spectrumsMatrix), dim(continuumsMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(nrow(spectrumsMatrix), length(sizesVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	if (is.null(windows)) 
	{
		windows <- loadDefaultFeWindows()
	}
	if (is.null(fitParameters)) 
	{
		fitParameters <- DEFAULT_FIT_PARAMETERS
	}
	return (cppFeFitRun(quasarcl, spectrumsMatrix, wavelengthsMatrix, errorsMatrix, 
							  continuumsMatrix, sizesVector, feTemplate, windows, 
							  fitParameters))
}



#' @export
rCalcFeTemplateMatrix <- function(quasarcl, wavelengthsMatrix, sizesVector, feTemplate,
								  fitParameters = NULL) 
{
	if(!isInitialized(quasarcl))
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(nrow(wavelengthsMatrix), length(sizesVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	if (is.null(fitParameters)) 
	{
		fitParameters <- DEFAULT_FIT_PARAMETERS
	}
	return (cppCalcFeTemplateMatrix(quasarcl, wavelengthsMatrix, sizesVector, feTemplate,
								    fitParameters))
}



#' @export
rCalcFeTemplateScaleRates <- function(quasarcl, spectrumsMatrix, feTemplateMatrix, sizesVector, 
									  fitParameters = NULL) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(dim(spectrumsMatrix), dim(feTemplateMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(nrow(spectrumsMatrix), length(sizesVector)))
	{
		stop("Invalid spectrums sizes vector length")
	}
	if (is.null(fitParameters)) 
	{
		fitParameters <- DEFAULT_FIT_PARAMETERS
	}
	return (cppCalcFeTemplateScaleRates(quasarcl, spectrumsMatrix, feTemplateMatrix, sizesVector, fitParameters))
}



#' @export
rCalcReducedChisqs <- function(quasarcl, fMatrix, yMatrix, errorsMatrix, sizesVector) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(dim(fMatrix), dim(yMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(dim(fMatrix), dim(errorsMatrix)))
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(nrow(fMatrix), length(sizesVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	return (cppCalcReducedChisqs(quasarcl, fMatrix, yMatrix, errorsMatrix, sizesVector))
}



#' @export
rCpuConvolve <- function(quasarcl, signalVector, kernelVector, same) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	return (cppCpuConvolve(quasarcl, signalVector, kernelVector, same))
}



#' @export
rReduceFeChisqs <- function(quasarcl, chisqsVector, sizesVector) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(length(chisqsVector), length(sizesVector))) 
	{
		stop("Vectors with different sizes")
	}
	return (cppReduceFeChisqs(quasarcl, chisqsVector, sizesVector))
}
