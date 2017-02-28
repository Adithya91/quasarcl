#' @useDynLib quasar
#' @importFrom Rcpp evalCpp



#' @export
rConvolve <- function(quasarcl, inputMatrix, sizesVector, filterVector) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(nrow(inputMatrix), length(sizesVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	if(!identical(nrow(inputMatrix), length(filterVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	return (cppConvolve(quasarcl, inputMatrix, sizesVector, filterVector))
}



#' @export
rCopyIfNotInf <- function(quasarcl, inputMatrix, filteredSize) 
{
	if(!isInitialized(quasarcl))
	{
		stop("QuasarCL is not initialized!")
	}
	return (cppCopyIfNotInf(quasarcl, inputMatrix, adjustSize(filteredSize)))
}



#' @export
rCountIfNotInf <- function(quasarcl, inputMatrix) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	return (cppCountIfNotInf(quasarcl, inputMatrix))
}



#' @export
rReglin <- function(quasarcl, xMatrix, yMatrix, sizesVector) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(dim(xMatrix), dim(yMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(nrow(xMatrix), length(sizesVector)))
	{
		stop("Invalid spectrums sizes vector length")
	}
	return (cppReglin(quasarcl, xMatrix, yMatrix, sizesVector))
}



#' @export
rReglinYax <- function(quasarcl, xMatrix, yMatrix, sizesVector) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(dim(xMatrix), dim(yMatrix)))
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(nrow(xMatrix), length(sizesVector)))
	{
		stop("Invalid spectrums sizes vector length")
	}
	return (cppReglinYax(quasarcl, xMatrix, yMatrix, sizesVector))
}



#' @export
rChisq <- function(quasarcl, fMatrix, yMatrix, errorsMatrix, sizesVector)
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
	return (cppChisq(quasarcl, fMatrix, yMatrix, errorsMatrix, sizesVector))
}



#' @export
rTrapz <- function(quasarcl, yMatrix, xMatrix, sizesVector)
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(dim(yMatrix), dim(xMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	if(!identical(nrow(yMatrix), length(sizesVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	return (cppTrapz(quasarcl, yMatrix, xMatrix, sizesVector))
}