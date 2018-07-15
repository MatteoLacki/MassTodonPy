require(Rcpp)
require(dplyr)
require(mzR)

combineSpectra <- function( path )
	path %>% openMSfile %>% peaks %>% lapply(
		function( spectrum )
		{
			spectrum 		<- as.data.frame( spectrum )
			names(spectrum) <- c('mz', 'intensity')
			spectrum %>% filter( intensity > 0 )
		}
	) %>%
	bind_rows %>%
	arrange(mz) 

combineSpectra.scanNo <- function( path )
{
	Data <-	path %>% openMSfile

	Map(
		function( spectrum, scanNo )
		{
			spectrum 		<- as.data.frame( spectrum )
			names(spectrum) <- c('mz', 'intensity')
			spectrum %>% 
			filter( intensity > 0 ) %>% 
			mutate( scanNo = scanNo )
		},
		Data %>% peaks,
		1:(Data%>%runInfo)$scanCount
	) %>%
	bind_rows %>%
	arrange(mz) 
}


sourceCpp('~/Dropbox/Science/MassSpectrometry/MassTodonTests/spectrumAveraging/DirkCluster.cpp')
DirkAverager <- function( spectrum, threshold = .008, diam = .1)
	spectrum %>%
	mutate( clustering = DirkCluster( mz, threshold, diam) ) %>%
	group_by( clustering ) %>%
	mutate(
		weights 	= intensity/sum(intensity),
		mzWeighted 	= mz * weights
	) %>%
	summarise(
		mz 			= sum(mzWeighted),
		intensity 	= sum(intensity),
		n 			= n()
	) %>% 
	select( -clustering ) %>%
	arrange( mz ) 

# require(Rcpp)
# require(dplyr)
sourceCpp('~/Dropbox/Science/MassSpectrometry/MassTodonTests/spectrumAveraging/HillCluster.cpp')
hillAverager <- function( spectrum, tol = .1 )
	spectrum %>%
		mutate(
		clustering = as.matrix(.) %>% HillCluster(tol)
	) %>%
	group_by( clustering ) %>%
	mutate(
		weights 	= intensity/sum(intensity),
		mzWeighted 	= mz * weights
	) %>%
	summarise(
		mz 			= sum(mzWeighted),
		intensity 	= sum(intensity),
		n 			= n()
	) %>% 
	select( -clustering ) %>%
	arrange( mz )


# data.frame(
# 	mz 			= 1:10,
# 	intensity 	= c(1,2,4,2,1,5,6,7,1,10)	 
# ) %>% hillAverager(., tol = 20 )

# dataPath 	<- '/Users/matteoLacki/Dropbox/MassTodon/Data/ubiquitin/July2015/12plus precursor/'
# setwd(dataPath)
