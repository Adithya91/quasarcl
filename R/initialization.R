#' @useDynLib quasar
#' @importFrom Rcpp evalCpp
	


#' @export
initialize <- function() 
{
	files <- dir(system.file("kernels", package = "quasar"), full.names = TRUE)
	if (identical(files, character(0))) 
	{
		stop("No OpenCL kernels found")
	}
	return (cppInitialize(files))
}



#' @export
isInitialized  <- function(quasarcl) 
{
	return (cppIsInitialized(quasarcl))
}