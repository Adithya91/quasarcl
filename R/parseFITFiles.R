#' @importFrom FITSio readFITS



parseFitFile <- function(file)
{
	fitFile <- readFITS(file)
	params <- parseHeader(fitFile)
	values <- parseData(fitFile)
	error <- parseError(fitFile)
	cat(".");
	return (list(params = params, values = values, error = error))
}



parseHeader <- function(fitFile)
{

	mjd <- as.integer(fitFile$hdr[which(fitFile$hdr == "MJD")+1])
	fiber <- as.integer(fitFile$hdr[which(fitFile$hdr == "FIBERID")+1])
	plate <- as.integer(fitFile$hdr[which(fitFile$hdr == "PLATEID")+1])
	b <- as.numeric(fitFile$hdr[which(fitFile$hdr == "COEFF0")+1])
	a <- as.numeric(fitFile$hdr[which(fitFile$hdr == "COEFF1")+1])
	z <- as.numeric(fitFile$hdr[which(fitFile$hdr == "Z")+1])
	ra <- as.numeric(fitFile$hdr[which(fitFile$hdr == "RAOBJ")+1])
	dec <- as.numeric(fitFile$hdr[which(fitFile$hdr == "DECOBJ")+1])
	type <- as.character(fitFile$hdr[which(fitFile$hdr == "OBJTYPE")+1])
	name <- coordinatesToName(ra, dec, "SDSS")

	return (list(mjd = mjd, fiber = fiber, plate = plate, b=b, a=a, z=z, ra=ra, dec=dec, type=type, name=name))
}



parseData <- function(fitFile)
{
	fluxo <- fitFile$imDat[,1]
	fluxo <- fluxo * 1.0e-17;
	return (fluxo)
}



parseError <- function(fitFile) 
{
	error <- fitFile$imDat[,3]
	error <- error * 1.0e-17;
	return (error)
}


