#' @useDynLib quasar
#' @importFrom Rcpp evalCpp



#' @export
rSimpleMAVG <- function(quasarcl, inputMatrix, windowWidth) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	return (cppSimpleMAVG(quasarcl, inputMatrix, windowWidth))
}



#' @export
rCenteredMAVG <- function(quasarcl, inputMatrix, sizesVector, windowWidth) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if(!identical(nrow(inputMatrix), length(sizesVector))) 
	{
		stop("Invalid spectrums sizes vector length")
	}
	return (cppCenteredMAVG(quasarcl, inputMatrix, sizesVector, windowWidth))
}

