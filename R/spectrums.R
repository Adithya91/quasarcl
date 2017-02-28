#' @useDynLib quasar
#' @importFrom Rcpp evalCpp
	


#' @export
rGenerateWavelengthsMatrix <- function(quasarcl, abzVector, size) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	return (cppGenerateWavelengthsMatrix(quasarcl, abzVector, size))
}



#' @export
rAddSpectrum <- function(quasarcl, wavelengthsMatrix, spectrumsMatrix, spectrumsSizesVector, 
						 toAddWavelengthsVector, toAddSpectrumVector) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(dim(spectrumsMatrix), dim(wavelengthsMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(nrow(wavelengthsMatrix), length(sizesVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	if(!identical(nrow(wavelengthsMatrix), length(toAddSpectrumVector))) 
	{
		stop("Invalid vector length")
	}
	if(!identical(nrow(wavelengthsMatrix), length(toAddWavelengthsVector))) 
	{
		stop("Invalid vector length")
	}
	return (cppAddSpectrum(quasarcl, wavelengthsMatrix, spectrumsMatrix, spectrumsSizesVector, 
						   toAddWavelengthsVector, toAddSpectrumVector))
}



#' @export
rFilterWithWavelengthWindows <- function(quasarcl, spectrumsMatrix, wavelengthsMatrix, errorsMatrix, 
										 sizesVector, windowsVector) 
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
	return (cppFilterWithWavelengthWindows(quasarcl, spectrumsMatrix, wavelengthsMatrix, 
										   errorsMatrix, sizesVector, windowsVector))
}



#' @export
rFilterNonpositive <- function(quasarcl, spectrumsMatrix, aMatrix, bMatrix, sizesVector) 
{
	if(!isInitialized(quasarcl))
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(dim(spectrumsMatrix), dim(aMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(dim(spectrumsMatrix), dim(bMatrix)))
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(nrow(spectrumsMatrix), length(sizesVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	return (cppFilterNonpositive(quasarcl, spectrumsMatrix, aMatrix, bMatrix, sizesVector))
}



#' @export
rFilterZeros <- function(quasarcl, spectrumsMatrix, aMatrix, bMatrix, sizesVector) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(dim(spectrumsMatrix), dim(aMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(dim(spectrumsMatrix), dim(bMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(nrow(spectrumsMatrix), length(sizesVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	return (cppFilterZeros(quasarcl, spectrumsMatrix, aMatrix, bMatrix, sizesVector))
}



#' @export
rFilterInfs <- function(quasarcl, spectrumsMatrix, aMatrix, bMatrix, sizesVector) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(dim(spectrumsMatrix), dim(aMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(dim(spectrumsMatrix), dim(bMatrix)))
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(nrow(spectrumsMatrix), length(sizesVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	return (cppFilterInfs(quasarcl, spectrumsMatrix, aMatrix, bMatrix, sizesVector))
}