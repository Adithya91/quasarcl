#' @export
drawGaussianRawData <- function(element, i, lambda) 
{
	if (element$fitParams[[i]][4] > 0) {
		lines(lambda, gaussian(lambda, element$fitParams[[i]][1] , element$fitParams[[i]][2], element$fitParams[[i]][3]), ylim=c(0,element$fitParams[[i]][1]))
	}
}



#' @export
drawGaussian <- function(element, lambda) 
{
	lines(lambda, gaussian(lambda, element$gaussianParams[1] , element$gaussianParams[2], element$gaussianParams[3], ylim=c(0,element$gaussianParams[1])))
}



#' @export
drawSpectrum <- function(wave, lambda, name) 
{
	plot(lambda,wave[1:length(lambda)],ylim=c(0, max(wave)*1.1), type='l', lwd=.2,lty=10,main=paste("Widmo kwazaru", name, sep=" "),xlab="dlugosc fali [A]",ylab="widmo [umowne jednostki]");
}


drawSpectrumWithPeaksRawData <- function(picturePath, spectrumsMatrix, wavelengthsMatrix, fitElements, qParams, sizes) {
	for (q in seq(1:length(sizes))) {
		jpeg(file=paste(picturePath, "widmo_", formatC(q, width=6, flag="0"),".jpg",sep=""),width = 500, height = 500, quality = 55, bg = "white")
		lambda <- wavelengthsMatrix[q,][1:sizes[q]]
		drawSpectrum1(spectrumsMatrix[q,][1:sizes[q]], lambda, qParams[[q]]$name)
		lapply(fitElements, drawGaussianRawData, q, lambda)
		dev.off();
	}
}
