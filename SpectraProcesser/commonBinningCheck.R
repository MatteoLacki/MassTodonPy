library(dplyr);library(ggplot2);library(tidyr);library(Rcpp);library(plotly)

# load spec
load("/Users/matteo/Documents/Science/MassTodon/MassTodonPy/spectrumPreprocessing/FRL_220715_ubi_952_ETD_40ms_01.Rdata")

# bitonicLocalMatrix
sourceCpp('~/Documents/Science/MassTodon/MassTodonPy/spectrumPreprocessing/bitonicSequences.cpp')

# is the binning similar? 

spec <- spec %>% arrange(mass) 
clusterInfo <- spec %>% mutate( clust = massDiffClustering(mass, thresh=2) ) %>% group_by(clust) %>% summarise(n=n())
clusterInfo$n %>% summary()

# are the estimated lower envel opes similar for 201 spectra?
massDiffMass <- function(spectrum){
	N <- nrow(spectrum)
	spectrum <- spectrum %>% arrange(mass)
	data.frame(
		mass 		= spectrum$mass[1:(N-1)],
		massDiff 	= spectrum$mass %>% diff,
		intensity 	= spectrum$intensity[1:(N-1)] 
	) %>% tbl_df
}

spec %>% split(.$specNo) %>% lapply(
	function(spec){
		clustering <- data.frame(massDiffMass(spec$mass, spec$intensity) )
		names(clustering) <- c('cluster', 'bigJump')	
	}
) 



# Here we compare the distribuant functions approximation.
probs <- seq(0,1,.01)

spec$mass %>% diff %>% quantile(probs=probs) %>%
data.frame(val = probs, quant=. ) %>% ggplot(
	aes(val,quant)
)+geom_line()+theme_minimal()+scale_y_log10()

D <-spec %>% split(.$specNo) %>% 
	lapply( function(d) 
		data.frame(
			prob 	= probs,
			quant 	= quantile(diff(d$mass), probs=probs)
		)
	) %>% bind_rows(.id='specNo') %>% mutate(specNo=as.integer(specNo)) %>% tbl_df

D %>% ggplot( aes(x=quant, y=prob, group=specNo, color=specNo) ) + geom_line() + theme_minimal() +
scale_x_log10()

# It seems that the curves are extremely similar, showing that there is something in the idea of 
# common thresholding. We can also compute the minimal distances within a spectrum, then join peaks,
# and finally 



MinPeakDistance <- spec %>% group_by( specNo ) %>% summarise( minDist = min( diff(sort(mass)) ) )
# all( MinPeakDistance$specNo > 0 ) # TRUE

minOfMins <- min(MinPeakDistance$minDist)

spec 		<- spec %>% mutate( clust = massDiffClustAbsolute(mass, minOfMins/10)) 
massSpreads <- spec %>% group_by(clust) %>% summarise( massSpread = max(mass)-min(mass), mass=mean(mass))



combinedSpec <- 
	spec %>% group_by(clust) %>% summarise( 
		n 			= n(), 	
		minMass 	= min(mass), 
		maxMass 	= max(mass), 
		meanMass 	= if(sum(intensity) > 0) crossprod(mass, intensity)/sum(intensity) else first(mass),
		intensity 	= sum(intensity) 
	)

combinedSpec %>% filter(meanMass>1712, meanMass<1716) %>%
ggplot(aes(x=minMass,y=intensity)) + geom_segment(aes(xend=meanMass,yend=0.0)) + geom_point() + theme_minimal()

combinedSpec %>% rename(mass=meanMass) %>% massDiffMass %>% select(mass, massDiff) %>%
filter( massDiff < 1.) %>% plot


# massSpreads %>% select(mass, massSpread) %>% filter(massSpread < 10^{-4}) %>% plot
spec %>% group_by(clust) %>% summarise( n=n(), intensity=mean(intensity)) %>% select(n,intensity) %>% 
group_by(n) %>% summarise(intensity=mean(intensity)) %>% plot(type='l')




