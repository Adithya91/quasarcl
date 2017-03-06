#' @useDynLib quasar
#' @importFrom Rcpp evalCpp




#' @title Logarytm dziesiętny
#' @description Oblicza logarytm dzięsiętny dla każdego elementu macierzy wejściowej
#' @param quasarcl (externalptr) Wskażnik na środowisko OpenCL
#' @param inputMatrix (matrix) Macierz wejściowa
#' @return (matrix) wynik
#' @export

rLog10 <- function(quasarcl, inputMatrix) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	return (cppLog10(quasarcl, inputMatrix))
}



#' @title Odejmowanie macierzy
#' @description Odejmuje od każdego elementu macierzy inputMatrix odpowiadający mu element
#' macierzy subtrahendMatrix
#' @param quasarcl (externalptr) Wskażnik na środowisko OpenCL
#' @param inputMatrix (matrix) Macierz wejściowa
#' @param subtrahendMatrix (matrix) Macierz odejmowana
#' @return (matrix) wynik 
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


#' @title Odejmowanie wartości skalarnej
#' @description Odejmuje od każdego elementu macierzy wejściowej skalar
#' @param quasarcl (externalptr) Wskażnik na środowisko OpenCL
#' @param inputMatrix (matrix) Macierz wejściowa
#' @param scalar (numeric) Skalar
#' @return (matrix) wynik 
#' @export

rMinusScalar <- function(quasarcl, inputMatrix, scalar) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	return (cppMinusScalar(quasarcl, inputMatrix, scalar))
}



#' @title Mnożenie kolumn macierzy przez wektor
#' @description Mnoży kolejne kolumny macierzy wejściowej przez kolejne wartości wektora 
#' @param quasarcl (externalptr) Wskażnik na środowisko OpenCL
#' @param inputMatrix (matrix) Macierz wejściowa
#' @param vector (numeric) Wektor
#' @return (matrix) wynik 
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



#' @title Transpozycja macierzy
#' @description Transponuje macierz wejściową 
#' @param quasarcl (externalptr) Wskażnik na środowisko OpenCL
#' @param inputMatrix (matrix) Macierz wejściowa
#' @return (matrix) wynik 
#' @export

rTranspose <- function(quasarcl, inputMatrix) 
{
	if(!isInitialized(quasarcl))
	{
		stop("QuasarCL is not initialized!")
	}
	return (cppTranspose(quasarcl, inputMatrix))
}



#' @title Dzielenie macierzy
#' @description Dzieli każdy element macierzy inputMatrix przez odpowiadający mu element 
#' macierzy divisorMatrix
#' @param quasarcl (externalptr) Wskażnik na środowisko OpenCL
#' @param inputMatrix (matrix) Macierz wejściowa
#' @param divisorMatrix (matrix) Macierz, przez której elementy dzielimy
#' @return (matrix) wynik 
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










