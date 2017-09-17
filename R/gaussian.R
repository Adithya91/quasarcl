#' @useDynLib quasar
#' @importFrom Rcpp evalCpp



#' @export
rParameterization <- function(quasarcl, quasars, options) 
{
	if(!isInitialized(quasarcl)) 
	{
		stop("QuasarCL is not initialized!")
	}
	
	spectrumsMatrix <- getSpectrumsMatrix(quasars)
	errorsMatrix <- getErrorsMatrix(quasars)
	sizesVector <- getSizesVector(quasars)
	abz <- getAbzParams(quasars)
	

	return (cppParameterization(quasarcl, spectrumsMatrix, errorsMatrix, sizesVector, abz, 
								options$spectralLines, options$continuumWindows, options$ampWavelength, 
								options$feWindows, options$feTemplate, options$fitParameters))
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
rFitGaussian <- function(quasarcl, yMatrix, xMatrix, sizesVector, resultsVector) 
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
		stop("Invalid spectrums sizes vector length")
	}
	if(!identical(nrow(yMatrix), length(resultsVector))) 
	{
		stop("Invalid results vector length")
	}
	return (cppFitGaussian(quasarcl, yMatrix, xMatrix, sizesVector, resultsVector))
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



#' @export
transformResults <- function(results, qParams, spectralLines) {
	transformedResults <- list();
	
	for (i in seq_along(qParams))
	{
		elementsFits <- mapply(getFitElement, results$elementsFits, spectralLines, i, SIMPLIFY=FALSE)
		elementsFits <- elementsFits[!is.na(elementsFits)]
		
		transformedResults[[i]] <- list(mjd = qParams[[i]]$mjd, 
										plate =  qParams[[i]]$plate, 
										fiber = qParams[[i]]$fiber, 
										continuumChisq = results$continuumChisqs[[i]], 
										continuumReglin = results$continuumReglin[[i]], 
										reglin = results$reglin[[i]], 
										feScaleRate = results$feScaleRates[[i]], 
										feWindowsSize = results$feWindowsSizes[[i]], 
										feWindowsReducedChisq = results$feWindowsReducedChisqs[[i]], 
										feFullReducedChisq = results$feFullReducedChisqs[[i]],
										feFullEW = results$feFullEWs[[i]], 
										feRangeReducedChisq = results$feRangeReducedChisqs[[i]], 
										feRangeEW = results$feRangeEWs[[i]], 
										elementsFits = elementsFits)
	}
	return(transformedResults)
}
