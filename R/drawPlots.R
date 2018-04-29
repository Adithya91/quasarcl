ASTRO_OBJ_SPEC_SIZE <- 4096

#' @export
drawGaussianRawData <- function(element, i, lambda) 
{
	if (element[4] > 0) {
		lines(lambda, gaussian(lambda, element[1] , element[2], element[3]), ylim=c(0,element[1]))
	}
}

#' @export
plotGaussian <- function(lambda, a, b, c) 
{
		plot(lambda, gaussian(lambda, a, b, c), ylim=c(0,a))
}


#' @export
drawGaussian <- function(element, lambda) 
{
	lines(lambda, gaussian(lambda, element$gaussianParams[1] , element$gaussianParams[2], element$gaussianParams[3], ylim=c(0,element$gaussianParams[1])))
}



#' @export
drawSpectrum <- function(wave, lambda, name) 
{
	plot(lambda,wave[1:length(lambda)],ylim=c(0, max(wave)*1.1), type='l', lwd=.2,lty=10,main=paste("Quasar ", name, sep=" "),xlab="wavelength [A]",ylab="flux (arbitrary units)");
}


#' @export
drawSpectrumWithPeaksRawData <- function(picturePath, spectrumsMatrix, wavelengthsMatrix, fitElements, qParams, sizes) {
	for (q in seq(1:length(sizes))) {
		jpeg(file=paste(picturePath, "spectrum_", formatC(q, width=6, flag="0"),".jpg",sep=""),width = 500, height = 500, quality = 55, bg = "white")
		lambda <- wavelengthsMatrix[q,][1:sizes[q]]
		drawSpectrum(spectrumsMatrix[q,][1:sizes[q]], lambda, qParams[[q]]$name)
		lapply(fitElements, drawGaussianRawData, q, lambda)
		dev.off();
	}
}

#' @export
drawChosenSpectrumWithPeaksRawData <- function(q, spectrumsMatrix, wavelengthsMatrix, fitElements, qParams, sizes) {
    lambda <- wavelengthsMatrix[q,][1:sizes[q]]
    drawSpectrum(spectrumsMatrix[q,][1:sizes[q]], lambda, qParams[[q]]$name)
    invisible(lapply(fitElements, drawGaussianRawData, q, lambda))
}

#' @export
rContinuumTest <- function(quasarcl, spectrumsMatrix, wavelengthsMatrix, errorsMatrix, sizes,
						 windows, ampWavelength) 
{
	filteredMatrices <- rFilterWithWavelengthWindows(quasarcl, spectrumsMatrix, wavelengthsMatrix, errorsMatrix, sizes, windows)
		
	filteredMatrices <- rFilterNonpositive(quasarcl, filteredMatrices$spectrumsMatrix, 
										   filteredMatrices$wavelengthsMatrix, filteredMatrices$errorsMatrix, sizes)

	newSizes <- rCountIfNotInf(quasarcl, filteredMatrices$spectrumsMatrix)
	maxSize <- max(newSizes)
	
	spectrumsMatrixFiltered <- rCopyIfNotInf(quasarcl, filteredMatrices$spectrumsMatrix, maxSize)
	wavelengthsMatrixFiltered <- rCopyIfNotInf(quasarcl, filteredMatrices$aMatrix, maxSize)
	errorsMatrixFiltered <- rCopyIfNotInf(quasarcl, filteredMatrices$bMatrix, maxSize)
	
	spectrumsMatrixFilteredCopy <- rLog10(quasarcl, spectrumsMatrixFiltered)
	wavelengthsMatrixFilteredCopy <- rLog10(quasarcl, wavelengthsMatrixFiltered)
	
	cReglinResults <- rReglin(quasarcl, wavelengthsMatrixFilteredCopy, spectrumsMatrixFilteredCopy, newSizes)

	if(ampWavelength > (.Machine$double.xmin ^ 10.001))
	{
		lampLog10 <- log10(ampWavelength)
		wavelengthsMatrixFilteredCopy <- rMinusScalar(quasarcl, wavelengthsMatrixFilteredCopy, lampLog10)
	}

	reglinResults <- rReglin(quasarcl, wavelengthsMatrixFilteredCopy, spectrumsMatrixFilteredCopy, newSizes)
	reglin <- rFixReglinResults(quasarcl, cReglinResults, reglinResults)
	
	continuumMatrixFiltered<- rCalcCw(quasarcl, wavelengthsMatrixFiltered, reglin$cReglinResults)
	chisqsFiltered <- rChisq(quasarcl, spectrumsMatrixFiltered, continuumMatrixFiltered, errorsMatrixFiltered,
							 newSizes)
	chisqsFiltered <- rReduceContinuumChisqs(quasarcl, chisqsFiltered, newSizes)

	cfunDcfun <- rCalcCfunDcfun(quasarcl, wavelengthsMatrix, reglin$cReglinResults, reglin$reglinResults)

	return(list(dcontinuumsMatrix = cfunDcfun$dcontinuums, 
				continuumsMatrix = cfunDcfun$continuums,
				continuumChisqs = chisqsFiltered,
				continuumReglin = reglin$reglinResults,
				reglin = reglin$cReglinResults))
}

#' @export
rFeFitTest <- function(quasarcl, spectrumsData, feTemplate, feWindows, fitParameters)
{
	spectrumsMatrix <- spectrumsData$spectrumsMatrix
	wavelengthsMatrix <- spectrumsData$wavelengthsMatrix
	errorsMatrix <- spectrumsData$errorsMatrix
	continuumsMatrix <- spectrumsData$continuumsMatrix
	sizes <- spectrumsData$sizes
	
	templateFeMatrix <- rCalcFeTemplateMatrix(quasarcl, wavelengthsMatrix, sizes, feTemplate, fitParameters)
	templateFeMatrixCopy <- templateFeMatrix

	if (!(fitParameters$isSubC))
	{
		spectrumsMatrix <- rMinusMatrix(quasarcl, spectrumsMatrix, continuumsMatrix)
	}

	if (fitParameters$fitType == "WIN" || fitParameters$fitType == "FWIN")
	{
		filteredMatrices <- rFilterWithWavelengthWindows(quasarcl, spectrumsMatrix, wavelengthsMatrix, errorsMatrix,
														 sizes, feWindows) 
		spectrumsMatrix <- filteredMatrices$spectrumsMatrix
		wavelengthsMatrix <- filteredMatrices$wavelengthsMatrix
		errorsMatrix <- filteredMatrices$errorsMatrix
		
		filteredMatrices <- rFilterInfs(quasarcl, spectrumsMatrix, continuumsMatrix, templateFeMatrix, sizes)
		spectrumsMatrix <- filteredMatrices$spectrumsMatrix
		continuumsMatrix <- filteredMatrices$aMatrix
		templateFeMatrix <- filteredMatrices$bMatrix
	}

	if (fitParameters$fitType == "FWIN")
	{
		filteredMatrices <- rFilterZeros(quasarcl, templateFeMatrix, wavelengthsMatrix, errorsMatrix, sizes)
		templateFeMatrix <- filteredMatrices$spectrumsMatrix
		wavelengthsMatrix <- filteredMatrices$aMatrix
		errorsMatrix <- filteredMatrices$bMatrix
		
		filteredMatrices <- rFilterInf(quasarcl, templateFeMatrix, continuumsMatrix, spectrumsMatrix, sizes)
		templateFeMatrix <- filteredMatrices$spectrumsMatrix
		continuumsMatrix <- filteredMatrices$aMatrix
		spectrumsMatrix <- filteredMatrices$bMatrix
	}
	
	spectrumsMatrixFiltered <- spectrumsMatrix;
	wavelengthsMatrixFiltered <- wavelengthsMatrix;
	errorsMatrixFiltered <- errorsMatrix;
	continuumsMatrixFiltered <- continuumsMatrix;
	templateFeMatrixFiltered <- templateFeMatrix;
	sizesFeWindows <- sizes;
	maxSpectrumSize <- ASTRO_OBJ_SPEC_SIZE
	
	if (fitParameters$fitType == "WIN" || fitParameters$fitType == "FWIN")
	{
		sizesFeWindows <- rCountIfNotInf(quasarcl, filteredMatrices$spectrumsMatrix)
		maxSize <- max(sizesFeWindows)
		
		spectrumsMatrixFiltered <- rCopyIfNotInf(quasarcl, spectrumsMatrix, maxSize)
		wavelengthsMatrixFiltered <- rCopyIfNotInf(quasarcl, wavelengthsMatrix, maxSize)
		errorsMatrixFiltered <- rCopyIfNotInf(quasarcl, errorsMatrix, maxSize)
		continuumsMatrixFiltered <- rCopyIfNotInf(quasarcl, continuumsMatrix, maxSize)
		templateFeMatrixFiltered <- rCopyIfNotInf(quasarcl, templateFeMatrix, maxSize)
	}

	scaleRates <- rCalcFeTemplateScaleRates(quasarcl, spectrumsMatrixFiltered, 
												templateFeMatrixFiltered, sizesFeWindows,
												fitParameters)
	
	templateFeMatrixFiltered <- rMultiplyCol(quasarcl, templateFeMatrixFiltered, scaleRates)
	templateFeMatrixCopy <- rMultiplyCol(quasarcl, templateFeMatrixCopy, scaleRates)
	
	reducedChisqsFiltered <- rCalcReducedChisqs(quasarcl, spectrumsMatrixFiltered, templateFeMatrixFiltered,
												errorsMatrixFiltered, sizesFeWindows)
												
												
	#pobranie pelnych danych dla dalszych obliczeń
	spectrumsMatrix <- spectrumsData$spectrumsMatrix
	wavelengthsMatrix <- spectrumsData$wavelengthsMatrix
	errorsMatrix <- spectrumsData$errorsMatrix
	continuumsMatrix <- spectrumsData$continuumsMatrix
	
	if (!(fitParameters$isSubC))
	{
		spectrumsMatrix <- rMinusMatrix(quasarcl, spectrumsMatrix, continuumsMatrix)
	}
	
	reducedChisqsFull <- rCalcReducedChisqs(quasarcl, spectrumsMatrix, templateFeMatrixCopy,
											errorsMatrix, sizes)
		
	temp <- rDivide(quasarcl, templateFeMatrixCopy, continuumsMatrix)
	ewsFull <- rTrapz(quasarcl, temp, wavelengthsMatrix, sizes)
	
	templateFeMatrix <- templateFeMatrixCopy
	
	if (fitParameters$fitType == "WIN" || fitParameters$fitType == "FWIN")
	{
		filteredMatrices <- rFilterWithWavelengthWindows(quasarcl, spectrumsMatrix, wavelengthsMatrix, 
														 errorsMatrix, sizes, list(fitParameters$feFitRange))

		spectrumsMatrix <- filteredMatrices$spectrumsMatrix
		wavelengthsMatrix <- filteredMatrices$wavelengthsMatrix
		errorsMatrix <- filteredMatrices$errorsMatrix 
		
		filteredMatrices <- rFilterInfs(quasarcl, spectrumsMatrix, continuumsMatrix, templateFeMatrix, 
										sizes)
										
		templateFeMatrix <- filteredMatrices$bMatrix
		continuumsMatrix <- filteredMatrices$aMatrix
		spectrumsMatrix <- filteredMatrices$spectrumsMatrix
		
		sizesFiltered <- rCountIfNotInf(quasarcl, spectrumsMatrix)
		maxSize <- max(sizesFiltered)
		
		spectrumsMatrix <- rCopyIfNotInf(quasarcl, spectrumsMatrix, maxSize)
		wavelengthsMatrix <- rCopyIfNotInf(quasarcl, wavelengthsMatrix, maxSize)
		errorsMatrix<- rCopyIfNotInf(quasarcl, errorsMatrix, maxSize)
		continuumsMatrix <- rCopyIfNotInf(quasarcl, continuumsMatrix, maxSize)
		templateFeMatrix <- rCopyIfNotInf(quasarcl, templateFeMatrix, maxSize)
	
		reducedChisqsFeRange <- rCalcReducedChisqs(quasarcl, spectrumsMatrix, templateFeMatrix,
												   errorsMatrix, sizesFiltered)
												   
		templateFeMatrix <- rDivide(quasarcl, templateFeMatrix, continuumsMatrix)
		ewsFeRange <- rTrapz(quasarcl, templateFeMatrix, wavelengthsMatrix, sizesFiltered)
		
	} else 
	{
		reducedChisqsFeRange <- reducedChisqsFull
		ewsFeRange <- ewsFull
	}
	
	return (list(feTemplateMatrix = templateFeMatrixCopy, 
				 feScaleRates = scaleRates, 
				 feWindowsSizes = sizesFeWindows, 
				 feWindowsReducedChisqs = reducedChisqsFiltered,
				 feFullReducedChisqs = reducedChisqsFull, 
				 feFullEWs = ewsFull, 
				 feRangeReducedChisqs = reducedChisqsFeRange, 
				 feRangeEWs = ewsFeRange))
}

#' @export
rFitElementTest <- function(quasarcl, spectrumsLinesMatrix, continuumsMatrix, wavelengthsMatrix, 
							errorsMatrix, sizesVector, element)
{
	spectrumsLinesMatrixCopy <- spectrumsLinesMatrix
	wavelengthsMatrixCopy <- wavelengthsMatrix
	sizesCopy <- sizesVector
	
	filteredMatrices <- rFilterWithWavelengthWindows(quasarcl, spectrumsLinesMatrix, wavelengthsMatrix, 
													 continuumsMatrix, sizesVector, list(element$range))
	
	sizes <- rCountIfNotInf(quasarcl, filteredMatrices$spectrumsMatrix)
	maxSize <- max(sizes)

	spectrumsLinesMatrix<- rCopyIfNotInf(ptr, filteredMatrices$spectrumsMatrix, maxSize)
	wavelengthsMatrix <- rCopyIfNotInf(ptr, filteredMatrices$wavelengthsMatrix, maxSize)
	continuumsMatrix <- rCopyIfNotInf(ptr, filteredMatrices$errorsMatrix, maxSize)
	
	fitGResults <- (rep(list(element$fitGuess), nrow(continuumsMatrix)))
	fitGResults <- rFitGaussian(quasarcl, spectrumsLinesMatrix, wavelengthsMatrix, sizes, fitGResults)
	
	gaussiansMatrix <- rCalcGaussian(quasarcl, wavelengthsMatrix, fitGResults, sizes)
	gaussiansMatrix <- rDivide(quasarcl, gaussiansMatrix, continuumsMatrix)
	ews <- rTrapz(quasarcl, gaussiansMatrix, wavelengthsMatrix, sizes)

	gaussianChisqs <- rCalcGaussianChisqs(quasarcl, wavelengthsMatrixCopy, spectrumsLinesMatrixCopy,
										 errorsMatrix, fitGResults, sizesCopy)
		
	gaussianFWHMs <- rCalcGaussianFWHM(quasarcl, fitGResults)
	
	return(list(fitParams = fitGResults,
				ews = ews,
				chisqs = gaussianChisqs,
				gaussianFWHMs = gaussianFWHMs))
	
}

#' @export
rParameterizationTest <- function(ptr, quasars, options) 
{
	#domyślne wartości - funkcje zbiorcze rContinuum, rFeFit oraz rParameterization ładują je we własnym zakresie,
	#nie trzeba tego robić ręcznie, jeżeli z nich korzystamy
	
	spectralLines <- options$spectralLines
	continuumWindows <- options$continuumWindows
	ampWavelength <- options$ampWavelength
	feWindows<- options$feWindows
	fitParameters <- options$fitParameters	
	feTemplate <- options$feTemplate

	errorsMatrix <- getErrorsMatrix(quasars)
	
	sizesVector <- getSizesVector(quasars)  

	spectrumsMatrix <- getSpectrumsMatrix(quasars)
	spectrumsMatrix <- rCenteredMAVG(ptr, spectrumsMatrix, sizesVector, 25)
	
	
	abz <- getAbzParams(quasars)
	wavelengthsMatrix <- rGenerateWavelengthsMatrix(ptr, abz, ASTRO_OBJ_SPEC_SIZE)
	wavelengthsMatrix <- rTranspose(ptr, wavelengthsMatrix)

	
	filteredMatrices <- rFilterZeros(ptr, errorsMatrix, spectrumsMatrix, wavelengthsMatrix, sizesVector)
	newSizes <- rCountIfNotInf(ptr, filteredMatrices$bMatrix)
	maxSize <- max(newSizes)

	spectrumsMatrix<- rCopyIfNotInf(ptr, filteredMatrices$aMatrix, ASTRO_OBJ_SPEC_SIZE)
	wavelengthsMatrix <- rCopyIfNotInf(ptr, filteredMatrices$bMatrix, ASTRO_OBJ_SPEC_SIZE)
	errorsMatrix <- rCopyIfNotInf(ptr, filteredMatrices$spectrumsMatrix, ASTRO_OBJ_SPEC_SIZE)

	spectrumsMatrixCopy <- spectrumsMatrix
	
	
	max_iter <- 3
	for (i in seq(1:max_iter)) 
	{
		continuumResults <- rContinuumTest(ptr, spectrumsMatrix, wavelengthsMatrix, errorsMatrix,
											newSizes, continuumWindows, ampWavelength)
		#continuumResults <- rContinuum(ptr, spectrumsMatrix, wavelengthsMatrix, errorsMatrix,
		#								newSizes, continuumWindows, ampWavelength)

		spectrumsMatrix <- spectrumsMatrixCopy
	
		spectrumsData <- list(spectrumsMatrix = spectrumsMatrix,
							  wavelengthsMatrix = wavelengthsMatrix,
							  errorsMatrix = errorsMatrix,
							  continuumsMatrix = continuumResults$continuumsMatrix,
							  sizes = newSizes)

		feResults <- rFeFitTest(ptr, spectrumsData, feTemplate, feWindows, fitParameters)
		#feResults <- rFeFit(ptr, spectrumsData, feTemplate, feWindows, fitParameters)

		if (i < max_iter)
		{
			spectrumsMatrix <- rMinusMatrix(ptr, spectrumsMatrix, feResults$feTemplateMatrix)
		}
	}

	spectrumsEmissionLines <- rMinusMatrix(ptr, spectrumsMatrix, feResults$feTemplateMatrix)
	spectrumsEmissionLines <- rMinusMatrix(ptr, spectrumsEmissionLines, continuumResults$continuumsMatrix)
	
	fitElementsResults <- lapply(spectralLines, rFitElementTest, quasarcl = ptr, 
																 spectrumsLinesMatrix = spectrumsEmissionLines, 
																 continuumsMatrix = continuumResults$continuumsMatrix, 
																 wavelengthsMatrix = wavelengthsMatrix, 
																 errorsMatrix = errorsMatrix, 
																 sizesVector = newSizes)
	
	return (list(continuumChisqs = continuumResults$continuumChisqs,
				 continuumReglin = continuumResults$continuumReglin,
				 reglin = continuumResults$reglin,
				 feScaleRates = feResults$feScaleRates,
				 feWindowsSizes = feResults$feWindowsSizes,
				 feWindowsReducedChisqs = feResults$feWindowsReducedChisqs,
				 feFullReducedChisqs = feResults$feFullReducedChisqs,
				 feFullEWs = feResults$feFullEWs,
				 feRangeReducedChisqs = feResults$feRangeReducedChisqs,
				 feRangeEWs = feResults$feRangeEWs,
				 elementsFits = fitElementsResults,
				 spectrumsEmissionLines = spectrumsEmissionLines,
				 fitElementsResults = fitElementsResults,
				 abz = abz,
				 newSizes = newSizes,
				 wavelengthsMatrix = wavelengthsMatrix,
				 spectrumsMatrix=spectrumsMatrix
				 ))

}

