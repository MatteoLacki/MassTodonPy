library(Rcpp)
library(mzR)

sourceCpp('~/Dropbox/Science/MassSpectrometry/MassTodonTests/spectrumAveraging/DirkAverager.cpp')

load('~/Dropbox/MassTodon/ProcessedData/ubiquitin/ETD_10:02:15.RData')

file <- MassTodons$files[1]
path <- 
'~/Dropbox/MassTodon/Data/ubiquitin/Dec_2014_orbitrap'


spectra <- openMSfile(
	paste0( path, '/', file )
)
	
spectra <- lapply(
	peaks( spectra ),
	function( spectrum ){

		spectrum <- as.data.frame( spectrum )
		names(spectrum) <- c('mz', 'intensity')

		return(
			spectrum %>% filter( intensity > 0 )
		)
	}
)

threshold 	<- .01
spectrum 	<- 	
	spectra %>% bind_rows %>% 
	arrange(
		mz
	) %>%
	mutate(
		clustering = clusterAllDirk( mz, threshold)
	) %>%
	group_by( clustering ) %>%
	summarise(
		mz 			= mean(mz),
		intensity 	= sum(intensity)
	) %>%
	select( mz, intensity ) 