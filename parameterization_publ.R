library(quasar)
library(latex2exp)
library(gdata)
library(reticulate)
ASTRO_OBJ_SPEC_SIZE <- 4096
'%&%' <- function(x, y)paste0(x,y)
mypath <- "/home/adithya/Astro/onefit/108"
qpath <- mypath
outpath <- "/home/adithya/Astro/astrout/"
dtpath <- "/usr/local/lib/R/site-library/quasar/data/"
#dir.create(mypath, showWarnings = TRUE, recursive = TRUE, mode = "0755")
setbasename<-"quasarset10_"
setparaname<-"qparamset10_"
N=10
i=1
quasarsetname<-setbasename %&% formatC(i, width=2, flag="0")
qparamsetname<-setparaname %&% formatC(i, width=2, flag="0")
quasardate<-"2018-04-30"
ptr <- initialize()
dbConnection <- getDbConnection()
start_time <- Sys.time()
quasars <- loadQuasarsFromFitFiles(qpath,1+(i-1)*N,N); 
length(quasars)
t1 <- Sys.time() - start_time; #t1
feTemplate <- loadIronEmissionTemplate(paste(dtpath,"iron_emission_temp",sep=""))
options <- list( spectralLines = loadDefaultSpectralLines(), 
                 continuumWindows = loadDefaultContinuumWindows(), 
                 ampWavelength =  loadDefaultAmpWavelength(), 
                 feWindows = loadDefaultFeWindows(),
                 fitParameters = loadDefaultFitParameters(),
                 feTemplate = feTemplate)
start_time <- Sys.time()
output <- rParameterization(ptr, quasars, options)
t4 <- Sys.time() - start_time; #t4
t4
#output
wavelengthsMatrix<-output$wavelengthsMatrix
#wavelengthsMatrix
spectrumsMatrix<-output$spectrumsMatrix
continuumsMatrix<-output$continuumsMatrix
feTemplatesMatrix<-output$feTemplatesMatrix
sizesVector<-output$sizes
outputfit<- transformResults(output, getParams(quasars), options$spectralLines)
#output<-NULL

chosen_q=1
drawSpectrumWithContinuum(chosen_q, quasars, wavelengthsMatrix, spectrumsMatrix, continuumsMatrix, sizesVector)

#chosen_q=1
drawSpectrumWithoutContinuumIron(chosen_q, quasars, wavelengthsMatrix, spectrumsMatrix, continuumsMatrix, feTemplatesMatrix, sizesVector)

#chosen_q=1
drawChosenPeaksComponents(chosen_q, quasars, outputfit, wavelengthsMatrix, spectrumsMatrix, continuumsMatrix, feTemplatesMatrix, sizesVector)

memorysort()

#plot(wavelengthsMatrix[q,][1:sizesVector[q]]);
q=1
plot(wavelengthsMatrix[q,][1:sizesVector[q]]);
