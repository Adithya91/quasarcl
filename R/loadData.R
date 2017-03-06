#' @importFrom utils file_test


#' @export
loadSpectralLines <- function(file)
{
	if(!file_test("-f", file))
	{
        stop("File with spectral lines does not exist")
    }
	df <- read.table(file, header = FALSE)
	colnames(df) <- c("name", "range0", "range1", "fitGuess0", "fitGuess1", "fitGuess2" )
	
	spectralLines <- apply(df[, 2:length(df)], 1, specLineToList)
	return (spectralLines)
}



#' @export
loadIronEmissionTemplate <- function(file) 
{
	if(!file_test("-f", file))
	{
        stop("File with template does not exist")
    }
	template <- read.table(file, header = FALSE)
	colnames(template) <- c("wavelengths", "values")

	return (list(wavelengths = template[, "wavelengths"], 
				 values = template[, "values"]))
}



#' @export
loadWindows <-function(file) 
{
	if(!file_test("-f", file))
	{
        stop("File with  windows does not exist")
    }
    df <- read.table(file, header = FALSE)
    return (split(as.matrix(df), seq(nrow(df))))
}



#' @export
loadQuasarsFromFitFiles <- function(path) 
{
	fitFiles <- list.files(path=path, pattern="\\.fit$", full.names = TRUE)
	if (identical(fitFiles, character(0))) 
	{
		stop("No .fit files found")
	}
	return (lapply(fitFiles, parseFitFile))
}



#' @export
getSpectrumsMatrix <- function(quasars) 
{
	spectrums <- lapply(lapply(quasars,  `[[`, 'values' ), 
						function(spectrum) { c(spectrum, rep(Inf, ASTRO_OBJ_SIZE - length(spectrum))) })
	return(do.call(rbind, spectrums))
}



#' @export
getErrorsMatrix <- function(quasars) 
{
	errors <- lapply(lapply(quasars,  `[[`, 'error' ), 
					 function(error) { c(error, rep(Inf, ASTRO_OBJ_SIZE - length(error))) })
	return(do.call(rbind, errors))
}



#' @export
getSizesVector <-function(quasars) 
{
	sizes <-sapply(lapply(quasars,  `[[`, 'values' ), length)
	return (sizes)
}



#' @export
getAbzParams <- function(quasars) 
{
	abz <- lapply(getParams(quasars), function(param) { c(param$a, param$b, param$z, 0)})
	return(abz)
}



getParams <- function(quasars) 
{
	params <-lapply(quasars, `[[`, 'params')
	return(params)
}



specLineToList <- function(elem) 
{
	return (list(range = c(elem[["range0"]], elem[["range1"]]), 
				 fitGuess = c(elem[["fitGuess0"]], elem[["fitGuess1"]], elem[["fitGuess2"]], 0)))
}