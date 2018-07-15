library(dplyr)
library(Rcpp)
library(ggplot2)
setwd('~/Dropbox/MassTodon/Data/ubiquitin/vary_WH_WV/test')

	# contains >>spectra<<
load('test.Rdata')


threshold 	<- .01
peaks 		<- spectra$peaks

S <- lapply(peaks, function(m){ 
		m <- as.data.frame(m)  
		colnames(m) <- c('mz', 'intensity') 

		return(m)
	} 
)

names(S) <- 1:length(S)

data <- 
	tidyr::unnest_( col = 'experimentNo', S) %>% 
	tbl_df %>%
	arrange( mz ) 

cppFunction('
	NumericVector clusterAllDirk( NumericVector masses, double tolerance)
	{
		int n = masses.size();
		NumericVector clustering(n);

		clustering[0] = 1;

		for( int i = 1; i < n; ++i){

			if( masses[i] - masses[i-1] >  tolerance ){ 
				clustering[i] = clustering[i-1] + 1;
			} else {
				clustering[i] = clustering[i-1];
			}
		}

		return clustering;
	}
')

data <-
	data %>% 
	arrange(
		mz
	) %>%
	bind_cols( 
		data.frame(clustering = clusterAllDirk(data$mz, threshold)) 
	) %>%  
	group_by( clustering ) %>%
	summarise(
		averageClusterMZ = mean(mz),
		totalClusterIntensity = sum(intensity) 
	) %>% 
	tbl_df

	# I take modulus 5 to use the same colours every 5 clusters.
averagedSpectraPlot <-
	data %>%
	mutate(
		clusteringMod5 = as.factor(clustering %% 5)
	) %>%
	filter(
		averageClusterMZ > 1000,
		averageClusterMZ < 2000
	) %>%
	ggplot(
		aes(
			x = averageClusterMZ,
			y = totalClusterIntensity,
			colour = clusteringMod5
		)
	) + 
	geom_segment(
		aes(
			xend = averageClusterMZ,
			yend = 0
		)
	) + 
	scale_colour_discrete(
		guide = FALSE
	) + 
	theme_bw()


# pdf( 
# 	file = 
# 		'~/Dropbox/MassTodon/Data/ubiquitin/vary_WH_WV/test.pdf',
# 	width 	= 200, 
# 	height 	= 10
# )
# 	averagedSpectraPlot
# dev.off()