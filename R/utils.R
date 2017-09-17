coordinatesToName <- function(ra, dec, prefix) 
{
	r <- c()
	d <- c()
	r[1] <- floor(ra / (360. / 24.))
	r[2] <- floor(((ra / (360. / 24.)) %% 1.) * 60.)
	r[3] = ((((ra / (360. / 24.)) %% 1.) * 60.) %% 1.) * 60.;
	
	d[1] = floor(abs(dec));
	d[2] = floor((abs(dec) %% 1.) * 60);
	d[3] = (((abs(dec) %% 1.) * 60) %% 1.) * 60.;

	sign = "+";
	if (dec < 0) {
		sign = "-";
	}

	sr <- c()
	sd <- c()
	for (i in 1:2) {
		sr[i] <- formatC(r[i], width=2, flag="0")
		sd[i] <- formatC(d[i], width=2, flag="0")
	}
	sr[3] <- formatC(r[3], width=5, flag="0", format="f", digits=2)
	sd[3] <- formatC(d[3], width=4, flag="0", format="f", digits=1)

	return(paste(prefix, paste(sr, collapse = ''), sign, paste(sd, collapse=''), sep=""))
}



adjustSize <- function(size) 
{
	newSize <- ifelse(size > 0, size, 64)
	remainder <- newSize %% 64
	if (remainder != 0)
	{
		newSize <- newSize + (64 - remainder)
	}
	if (newSize < size)
	{
		stop("Invalid spectrum size")
	}
	return (newSize)
}



getNObj <- function(someList, N, i) 
{
	ifelse(((i+N)<=length(someList)), last<- i+N-1, last <- length(someList))
	return (someList[i:last])
}




specLineToList <- function(elem) 
{
	return (list(name=elem[["name"]], range = c(as.numeric(elem[["range0"]]), as.numeric(elem[["range1"]])), 
				 fitGuess = c(as.numeric(elem[["fitGuess0"]]), as.numeric(elem[["fitGuess1"]]), as.numeric(elem[["fitGuess2"]]), 0)))
}



downloadFitFile <- function(quasarInfo, path)
{
	mjd <-formatC(as.integer(quasarInfo[["V47"]]), width=5, flag="0")
	plate <- formatC(as.integer(quasarInfo[["V48"]]), width=4, flag="0")
	fiber <- formatC(as.integer(quasarInfo[["V49"]]), width=3, flag="0")
	
	filename <- paste(paste("spSpec", mjd, plate, fiber, sep="-"), "fit", sep=".")
	address <- paste(SERVER_ADDRESS, plate, "1d", filename, sep="/")
	file <- paste(path, filename, sep="")
	download.file(url=address, destfile=file, method="curl")
}



gaussian<-function(x, a, b, c) {
	return (a * exp((-0.5) * ((x - b)*(x - b)) / (c*c)))
}



getFitElement <- function(element, line, q)
{
	if (element$fitParams[[q]][4] > 0)
	{
		return(list(name = line$name, 
					gaussianParams = c(element$fitParams[[q]][1], element$fitParams[[q]][2], element$fitParams[[q]][3]), 
					gaussianFWHM = element$gaussianFWHMs[[q]], 
					chisq = element$chisqs[[q]],
					ew = element$ews[[q]]))
	}
	else return(NA)
}
