#' @useDynLib quasar
#' @importFrom Rcpp evalCpp



#' @export
rContinuum <- function(quasarcl, spectrumsMatrix, wavelengthsMatrix, errorsMatrix, sizesVector, 
					   windows = NULL, ampWavelength = NULL) 
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
	if(!identical(nrow(spectrumsMatrix), length(sizesVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	if (is.null(windows)) 
	{
		windows <- loadDefaultContinuumWindows()
	}
	if (is.null(ampWavelength)) 
	{
		ampWavelength <- DEFAULT_AMP_WAVELENGTH
	}
	return (cppContinuum(quasarcl, spectrumsMatrix, wavelengthsMatrix, errorsMatrix, sizesVector, 
						 windows, ampWavelength))
}



#' @export
rFixReglinResults <- function(quasarcl, cReglinResultsVector, reglinResultsVector)
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(length(reglinResultsVector), length(cReglinResultsVector))) 
	{
		stop("Vectors with different sizes")
	}
	return (cppFixReglinResults(quasarcl, cReglinResultsVector, reglinResultsVector))
}



#' @export
rCalcCfunDcfun <- function(quasarcl, wavelengthsMatrix, cReglinResultsVector, reglinResultsVector) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(nrow(wavelengthsMatrix), length(cReglinResultsVector))) 
	{
		stop("Invalid continuum reglin results vector length")
	}
	if(!identical(nrow(wavelengthsMatrix), length(reglinResultsVector))) 
	{
		stop("Invalid dcontinuum reglin results vector length")
	}
	return (cppCalcCfunDcfun(quasarcl, wavelengthsMatrix, cReglinResultsVector, reglinResultsVector))
}



#' @export
rCalcCw <- function(quasarcl, wavelengthsMatrix, cReglinResultsVector) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(nrow(wavelengthsMatrix), length(cReglinResultsVector))) 
	{
		stop("Invalid continuum reglin results vector length")
	}
	return (cppCalcCw(quasarcl, wavelengthsMatrix, cReglinResultsVector))
}



#' @export
rReduceContinuumChisqs <- function(quasarcl, chisqsVector, sizesVector) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(length(chisqsVector), length(sizesVector))) 
	{
		stop("Vectors with different sizes")
	}
	return (cppReduceContinuumChisqs(quasarcl, chisqsVector, sizesVector))
}
