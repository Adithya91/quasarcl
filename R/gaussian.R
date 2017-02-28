#' @useDynLib quasar
#' @importFrom Rcpp evalCpp



#' @export
rParameterization <- function(quasarcl, quasars, feTemplate, 
							  spectralLines = NULL, continuumWindows = NULL, 
							  ampWavelength = NULL, feWindows = NULL, fitParameters = NULL) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	
	spectrumsMatrix <- getSpectrumsMatrix(quasars)
	errorsMatrix <- getErrorsMatrix(quasars)
	sizesVector <- getSizesVector(quasars)
	abz <- getAbzParams(quasars)
	
	if (is.null(spectralLines)) 
	{
		spectralLines <- loadDefaultSpectralLines()
	}
	if (is.null(continuumWindows)) 
	{
		continuumWindows <- loadDefaultContinuumWindows()
	}
	if (is.null(ampWavelength)) 
	{
		ampWavelength <- DEFAULT_AMP_WAVELENGTH
	}
	if (is.null(feWindows)) 
	{
		feWindows<- loadDefaultFeWindows()
	}
	if (is.null(fitParameters)) 
	{
		fitParameters <- DEFAULT_FIT_PARAMETERS
	}

	return (cppParameterization(quasarcl, spectrumsMatrix, errorsMatrix, sizesVector, abz, 
								spectralLines, continuumWindows, ampWavelength, 
								feWindows, feTemplate, fitParameters))
}



#' @export
rFitElement <- function(quasarcl, specLinesMatrix, continuumsMatrix, wavelengthsMatrix, errorsMatrix, 
						sizesVector, element)
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(dim(specLinesMatrix), dim(continuumsMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(dim(specLinesMatrix), dim(wavelengthsMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(dim(specLinesMatrix), dim(errorsMatrix)))
	{
		stop("Matrices with different dimensions")
	}	
	if(!identical(nrow(specLinesMatrix), length(sizesVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	
	return (cppFitElement(quasarcl, specLinesMatrix, continuumsMatrix, wavelengthsMatrix, errorsMatrix, 
						  sizesVector, element))
}



#' @export
rFitGaussian <- function(quasarcl, yMatrix, xMatrix, sizesVector) 
{
	if(!isInitialized(quasarcl))
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(dim(xMatrix), dim(yMatrix)))
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(nrow(yMatrix), length(sizesVector))) 
	{
		stop("Matrices with different dimensions")
	}
	return (cppFitGaussian(quasarcl, yMatrix, xMatrix, sizesVector))
}



#' @export
rCalcGaussian <- function(quasarcl, xMatrix, gaussianParamsVector, sizesVector) 
{
	if(!isInitialized(quasarcl))
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(nrow(xMatrix), length(gaussianParamsVector)))
	{
		stop("Invalid spectrums sizes vector length")
	}
	if(!identical(nrow(xMatrix), length(sizesVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	return (cppCalcGaussian(quasarcl, xMatrix, gaussianParamsVector, sizesVector))
}



#' @export
rCalcGaussianChisqs <- function(quasarcl, xMatrix, yMatrix, errorsMatrix, gaussianParamsVector, 
								sizesVector) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(dim(xMatrix), dim(yMatrix)))
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(dim(xMatrix), dim(errorsMatrix)))
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(nrow(xMatrix), length(gaussianParamsVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	if(!identical(nrow(xMatrix), length(sizesVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	return (cppCalcGaussianChisqs(quasarcl, xMatrix, yMatrix, errorsMatrix, gaussianParamsVector, 
								  sizesVector))
}



#' @export
rCalcGaussianFWHM <- function(quasarcl, gaussianParamsVector) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	return (cppCalcGaussianFWHM(quasarcl, gaussianParamsVector))
}

