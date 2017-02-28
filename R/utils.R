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
	newSize = size
	remainder = newSize %% 64
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

