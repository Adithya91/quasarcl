DEFAULT_AMP_WAVELENGTH <- 3000.0

DEFAULT_FIT_PARAMETERS <- list(fwhmn = 1600.0, fwhmt = 900.0, feScaleRate = 1.0, 
							   feFitRange = c(2200.0, 2650.0), isSubC = FALSE, fitType="WIN")

ASTRO_OBJ_SIZE <- 4096


#' @export   
loadDefaultContinuumWindows <- function() 
{
	file <- system.file("data", "cont_windows",  package = "quasar")
	return (loadWindows(file))
}



#' @export
loadDefaultFeWindows <- function() 
{
	file <- system.file("data", "iron_emission_windows",  package = "quasar")
	return (loadWindows(file))
}



#' @export
loadDefaultSpectralLines <-function() 
{
	file <- system.file("data", "spectral_lines",  package = "quasar")
	return (loadSpectralLines(file))
}



#' @export
loadDefaultFitParameters <- function()
{
	return (DEFAULT_FIT_PARAMETERS)
}



#' @export
loadDefaultAmpWavelength <- function()
{
	return (DEFAULT_AMP_WAVELENGTH)
}
