#' @useDynLib quasar
#' @importFrom Rcpp evalCpp



#' @export
rLog10 <- function(quasarcl, inputMatrix) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	return (cppLog10(quasarcl, inputMatrix))
}



#' @export
rMinusMatrix <- function(quasarcl, inputMatrix, subtrahendMatrix) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if (!identical(dim(inputMatrix), dim(subtrahendMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	return(cppMinusMatrix(quasarcl, inputMatrix, subtrahendMatrix))
}



#' @export
rMinusScalar <- function(quasarcl, inputMatrix, scalar) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	return (cppMinusScalar(quasarcl, inputMatrix, scalar))
}



#' @export
rMultiplyCol <- function(quasarcl, inputMatrix, vector) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if (!identical(nrow(inputMatrix), length(vector))) 
	{
		stop("Invalid vector size")
	}
	return (cppMultiplyCol(quasarcl, inputMatrix, vector))
}



#' @export
rTranspose <- function(quasarcl, inputMatrix) 
{
	if(!isInitialized(quasarcl))
	{
		stop("QuasarCL is not initialized!")
	}
	return (cppTranspose(quasarcl, inputMatrix))
}



#' @export
rDivide <- function(quasarcl, inputMatrix, divisorMatrix) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	if (!identical(dim(inputMatrix), dim(divisorMatrix))) 
	{
		stop("Matrices with different dimensions")
	}
	return (cppDivide(quasarcl, inputMatrix, divisorMatrix))
}










