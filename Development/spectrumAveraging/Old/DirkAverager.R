require(Rcpp)
require(dplyr)
require(mzR)

# cppFunction('
# NumericVector clusterAllDirk(
#     NumericVector masses,
# 	double tolerance
# ){
# 	int n = masses.size();
# 	NumericVector clustering(n);

# 	clustering[0] = 1;

# 	for( int i = 1; i < n; ++i){

# 		if( masses[i] - masses[i-1] >  tolerance ){
# 			clustering[i] = clustering[i-1] + 1;
# 		} else {
# 			clustering[i] = clustering[i-1];
# 		}
# 	}

# 	return clustering;
# }
# ')

cppFunction('
NumericVector clusterAllDirk(
    NumericVector masses,
	double tolerance,
	double diam
){
	int n = masses.size();
	NumericVector clustering(n);

	clustering[0] = 1;
	double currDiam = 0.0;

	for( int i = 1; i < n; ++i){
		if( 
			(masses[i] - masses[i-1] >  tolerance) ||
			( currDiam > diam )
		){
			clustering[i] = clustering[i-1] + 1;
			currDiam = 0.0;
		} else {
			clustering[i] = clustering[i-1];
			currDiam += masses[i] - masses[i-1];
		}
	}

	return clustering;
}
')

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


DirkAverager <- function( spectrum, threshold = .008, diam = .1)
	spectrum %>%
	mutate( clustering = clusterAllDirk( mz, threshold, diam) ) %>%
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
	arrange( mz ) 


# dataPath 	<- '/Users/matteoLacki/Dropbox/MassTodon/Data/ubiquitin/July2015/12plus precursor/'
# setwd(dataPath)
