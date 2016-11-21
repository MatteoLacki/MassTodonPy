library(readMzXmlData)
library(dplyr)
library(scales)
library(ggplot2)
library(plotly)
library(tidyr)
library(Rcpp)
library(grDevices) # chull

# install.packages('plotly')
############################################################

plotSpec <- function(data){
	specPlot <- data %>% 
	ggplot(aes(x=mass,y=intensity))+
	geom_segment(aes(xend=mass,yend=0.0))
	return(specPlot)
}

almostLogTransform <- function(x) ifelse(x > 0, 1+log(x), x)

############################################################

Spectrum <- readMzXmlData::readMzXmlFile(
	mzXmlFile 		= "~/Documents/Science/MassTodon/MassTodonPy/data/FRL_220715_ubi_952_ETD_40ms_01.mzXML",
	removeMetaData 	= FALSE,
	verbose 		= TRUE
) 

######## Number of peaks above a given intensity threshold.
threshold_on_intensity <- -1
Spectrum %>% sapply(function(x) sum(x$spectrum$intensity > threshold_on_intensity) )


######## Spectra of the World - unite!
spec <- lapply(
	Spectrum,
	function( info ) as.data.frame( info$spectrum ) %>% arrange(mass)
) %>% bind_rows(.id='specNo') 		%>% 
filter(intensity > threshold_on_intensity) %>% 
mutate(specNo=as.integer(specNo)) 	%>% tbl_df 
# save(spec, file='FRL_220715_ubi_952_ETD_40ms_01.Rdata')


metaDatas <- Spectrum %>% lapply(
	function(x)
		x$metaData[8:27] %>% data.frame
) %>% bind_rows %>% tbl_df

# metaDatas$precursor.precursorCharge

# metaDatas %>% data.frame %>% head
# metaDatas %>% select( retentionTime, totIonCurrent) %>% 
# plot(type='l')

############################################################

spec %>% group_by(mass) %>% summarise(n=n()) %>% arrange(desc(n)) %>% select(n) %>% unlist %>% summary

roundedSpec <- spec %>%
	mutate( mass=round(mass,2) ) %>% group_by( mass ) %>%
	summarise( intensity=sum(intensity) )

# A <-spec %>% group_by(specNo) %>% summarise(n=n()) %>% 
# 	mutate( specNo=as.integer(specNo) ) 
# summary(A$n)

SpecVariability <- 
	spec %>% filter(mass>1712,mass<1716) %>% 
	mutate(mass = round(mass,2)) %>%
	group_by(mass, specNo) %>%
	summarise(intensity = sum(intensity) )%>%
	mutate( specNo=as.numeric(specNo)/201  )%>%
	ggplot( aes( x=mass, y=specNo, alpha=intensity ))+
	geom_point(size=.1)+theme_minimal()
# ggsave('SpecVariability2.pdf')

spec %>% filter(mass>1712,mass<1716) %>%  
mutate(
	mass 		= mass+rnorm(n(), sd=.005), 
	intensity 	= intensity+rnorm(n(),sd=.005)
) %>%
ggplot(aes( x=mass, y=intensity, color=specNo ))+
geom_point(size=.1)

spec %>% plotSpec()
A <- spec %>% filter(specNo==100, mass>1713.5, mass<1716) %>% 
plotSpec()
A+theme_bw()+geom_point()
ggplotly()

#### Isn't that all a gaussian curve evaluated?
spec %>% plotSpec()
A <- spec %>% filter(specNo==100, mass>1714.43, mass<1714.64) %>% 
plotSpec()
A+theme_bw()+geom_point()
ggplotly()

D 		<- spec %>% filter(specNo==100, mass>1714.43, mass<1714.64)
mStart 	<- crossprod( D$mass, D$intensity/sum(D$intensity))
dStart 	<-(max(D$mass)-min(D$mass))/4.5
maxPeakInd 	<- which.max(D$intensity)
maxPeakM	<- D$mass[maxPeakInd]
maxPeakI 	<- D$intensity[maxPeakInd]
nStart 		<- exp(log(maxPeakI) - dnorm(maxPeakM, mean=mStart, sd=dStart ,log=TRUE))
# fitting <- function(x){
# 	m <- x[1]; d <- x[2]; n <- x[3]
# 	crossprod( dnorm( D$mass, mean=m, sd=d )*n - D$intensity)	
# }

E <- D
bind_rows(
	D %>% mutate(tag = 'experiment'),
	data.frame(
		mass 		= D$mass,
		intensity 	= dnorm( D$mass, mean=mStart, sd=dStart )*nStart,
		tag	 		= 'theory'
	)
) %>% 
ggplot(aes(mass,intensity,color=tag,group=tag))+
geom_line()+
geom_point()+
theme_minimal()+
ggtitle('There must be some local Gaussian windowing.')
# ggsave('localGaussianWindowing.pdf')

#################################################################

# Checking if smallest peak distances are functions of m/z. This can be done so much brighter...
P<- spec %>% filter(
		mass>2855, 
		mass<2859,
		specNo==100
	) %>% 
	mutate(specReal=as.numeric(specNo)/201) %>%
	ggplot(aes( mass, intensity,color=specReal,group=specNo))+
	geom_line()+
	geom_point()+
	theme_minimal() + scale_colour_gradient(low = "white", high = "black")
# ggsave('art_spec.png')

summarise2 <- function(D) D %>% select(mass) %>% unlist %>% diff %>% summary
spec %>% filter(mass>942,  mass<944, 	specNo==100) %>% summarise2
spec %>% filter(mass>2131, mass<2131.5, specNo==100) %>% summarise2
spec %>% filter(mass>2855, mass<2859, 	specNo==100) %>% summarise2
ggplotly()

massDiffMass <- function(spectrum){
	N <- nrow(spectrum)
	spectrum <- spectrum %>% arrange(mass)
	data.frame(
		mass 		= spectrum$mass[1:(N-1)],
		massDiff 	= spectrum$mass %>% diff,
		intensity 	= spectrum$intensity[1:(N-1)] 
	) %>% tbl_df
}

plotSpec2 <- function(data) data %>% ggplot(aes(x=mass, y=massDiff))+geom_line()+theme_minimal()
plotSpec3 <- function(data) data %>% ggplot(aes(x=mass, y=massDiff))+geom_point()+theme_minimal()

### HA: the effect cannot be visually spotted if we use to big masses.
spec100massVSmassDiff <- spec %>% filter(specNo==100) %>% massDiffMass 
spec100massVSmassDiff %>% filter( massDiff < .1 ) %>% plotSpec3
spec100massVSmassDiff %>% filter( massDiff < .5 ) %>% plotSpec3
spec100massVSmassDiff %>% filter( massDiff < 1.2 )%>% plotSpec3

spec %>% split(.$specNo) %>% lapply(massDiffMass) %>% bind_rows(.id='specNo') %>%
filter( massDiff < 1) %>% ggplot(aes(x=mass, y=massDiff, color=as.numeric(specNo)))+geom_point()+theme_minimal()
ggplotly()

minDiffs %>% 


# How many peaks are there that are farther 
# than .02 than their right neighbor?

x 		<- rnorm(1000)
neighNo <- function(x, thr=.02){
	diffs <- sort(x) %>% diff 
	sum(diffs>thr)	
}
spec %>% group_by(specNo) %>% summarise(N=neighNo(mass), n=n()) %>%
gather("number","val",2:3)%>%
ggplot(aes(specNo,val,color=number))+geom_line()

spec %>% group_by(specNo) %>% summarise(N=neighNo(mass,thr=.2), n=n()) %>%
mutate(ratio = N/n) %>% select(specNo, ratio) %>% plot(type='l')


minDistNo <- function(spectrum,relative=FALSE){
	diffs <- hist(
		x 		= diff(sort(spectrum$mass)),
		breaks 	= c(seq(0,2,length.out=1000),Inf), 
		plot 	= FALSE
	)
	N <- sum(diffs$counts)
	y <- N - cumsum(diffs$counts)
	if(relative) y <- y/N
	y <- y[-length(y)]
	x <- diffs$mids
	x <- x[-length(x)]
	data.frame(massDiff=x,cnt=y)
}

massDiffs <- spec %>% split(.$specNo)%>% lapply(minDistNo) %>%
	bind_rows(.id='specNo') %>% tbl_df() %>% 
	mutate(specNo=as.integer(specNo))
massDiffs2 <- spec %>% split(.$specNo)%>% lapply(minDistNo,relative=TRUE) %>%
	bind_rows(.id='specNo') %>% tbl_df() %>% 
	mutate(specNo=as.integer(specNo))

massDiffs %>% 
filter(massDiff<0.082)%>%
ggplot(aes(massDiff,cnt,color=specNo/201,group=specNo))+geom_line()+
scale_x_continuous(breaks=seq(0,2,by=.01))+
theme_minimal()+
theme(axis.text.x = element_text(angle = 90, hjust = 1))+
ggtitle('Number of clusters as depending on the minimal distance between peaks in cluster.')+
xlab('Minimal Distance between Peaks in a Cluster')+
ylab('Number of Clusters')
# ggsave('clustersNoVSminimalPeakDistance.pdf')


massDiffs2 %>% 
filter(massDiff<0.082)%>%
ggplot(aes(massDiff,cnt,color=specNo/201,group=specNo))+geom_line()+
scale_x_continuous(breaks=seq(0,2,by=.01))+
scale_y_continuous(breaks=seq(0,1,by=.02))+
theme_minimal()+
theme(axis.text.x = element_text(angle = 90, hjust = 1))+
ggtitle('Number of clusters to Number of Peaks as depending on the minimal distance between peaks in cluster.')+
xlab('Minimal Distance between Peaks in a Cluster')+
ylab('Number of Clusters to Number of Peaks')
ggsave('clustersNo2totalPeaksNoVSminimalPeakDistance.pdf')

######################################################



SpecVariability
ggplotly()

forPlot %>% mutate( y= as.numeric(specNo)/201) %>%
ggplot(aes(x=M, y=y, alpha=f(totInt)))+geom_point(size=.1)+theme_minimal()

forPlot %>% mutate( y= as.numeric(specNo)/201) %>% 
filter( M > 1000, M < 1100) %>%
ggplot(aes(x=M, y=y, alpha=f(totInt)))+geom_point(size=.1)+theme_minimal()



spec %>% filter(specNo == 2) %>% select(mass, intensity) %>% 
filter(mass > 1500, mass < 1800) %>%
plot(type='h')


spec %>% filter(specNo == 2) %>% 
ggplot(aes(x=mass,y=intensity))+
geom_segment(aes(xend=mass,yend=0.0))+theme_minimal()


ggplotly()

############################################################
# Naive bitonic clustering.

betterData 	<- spec %>% split(.$specNo) %>% lapply(massDiffMass) %>% bind_rows(.id='specNo') 
betterData %>% filter( massDiff < 1.) %>%
ggplot(aes(x=mass,y=massDiff,color=as.numeric(specNo))) + theme_minimal() + geom_point()

betterData %>% filter( massDiff < .25) %>%
ggplot(aes(x=mass,y=massDiff,color=as.numeric(specNo))) + theme_minimal() + geom_point()


betterData %>% filter( massDiff < 10, massDiff > 1) %>%
ggplot(aes(x=mass,y=massDiff,color=as.numeric(specNo))) + theme_minimal() + geom_point()


betterData %>% filter( specNo==100, massDiff < 1) %>%
ggplot(aes(x=mass,y=massDiff,color=as.numeric(specNo))) + theme_minimal() + geom_point()


betterData


testSpec <- spec %>% filter(specNo==100)
sourceCpp('~/Documents/Science/MassTodon/MassTodonPy/spectrumPreprocessing/bitonicSequences.cpp')

testSpec <- testSpec %>% arrange(mass) %>% mutate( cluster = bitonicLocal(mass, intensity) ) 

testSpec %>% filter(mass>1712,mass<1716) %>% mutate( cluster = factor(cluster) ) %>%
ggplot(aes(x=mass,y=intensity,color=cluster)) + geom_segment(aes(xend=mass,yend=0.0)) + theme_minimal() + geom_point()
# ggplotly()

# testSpec$intensity %>% summary
# testSpec %>% filter( intensity > 7934 ) %>% mutate( cluster = factor(cluster) ) %>%
# ggplot(aes(x=mass,y=intensity,color=cluster)) + geom_segment(aes(xend=mass,yend=0.0)) + theme_minimal() + geom_point() +  theme(legend.position="none")



massDiffMass2 <- function(Dupa){
	N <- nrow(Dupa)
	data.frame(
		mass 		= Dupa$mass[1:(N-1)],
		massDiff 	= diff(Dupa$mass),
		intensity 	= Dupa$intensity[1:(N-1)],
		jump  		= Dupa$jump[1:(N-1)],
		clustChange = diff(Dupa$cluster),
		cluster 	= Dupa$cluster[1:(N-1)]
	) %>% tbl_df	
}

bigJumps <- bitonicLocal2(testSpec$mass, testSpec$intensity) 

data_clustered <- testSpec %>% select(-cluster) %>%
	bind_cols(
		bigJumps %>% data.frame %>% tbl_df %>% rename( cluster=X1, bigJump=X2 )
	) %>% 
	mutate( jump = factor(ifelse(bigJump==0,'small','big')) ) %>%
	massDiffMass2

data_clustered2 <- data_clustered %>% filter(massDiff < 1.) %>% 
	mutate( totChange = interaction(jump, ifelse(clustChange==0,'sameCluster','changedCluster')) )

data_clustered2 %>% ggplot(
	aes( x=mass, y=massDiff, color=totChange,
		 size=ifelse(as.character(jump)=='small',.2,1)) )+
	geom_point(alpha=.5)+theme_dark()
# ggplotly()
data_clustered2 %>% group_by( jump ) %>% summary()

# Only the envelope points are very close. 
data_clustered2 %>% ggplot(
	aes( x=mass, y=massDiff, color=totChange,
		 size=ifelse(as.character(totChange)=='small.changedCluster',.2,1)) )+
	geom_point(alpha=.5)+theme_dark()


data_clustered2 %>% filter( jump=='small' ) %>% 
ggplot(
	aes( x=mass, y=massDiff, size=ifelse(clustChange==0,.3,2))
)+
geom_point(alpha=.5)+theme_minimal()

data_clustered %>% filter(massDiff < 1.) %>% summary()

# Strange that the big change does not happen together with big jump. Look into it.
data_clustered %>% filter(massDiff < 1.) %>% mutate(totChange=interaction(jump, ifelse(clustChange==0,'sameCluster','changedCluster'))) %>% group_by(totChange) %>% summarise(n=n())




betterData %>% filter( specNo==100, mass >310, mass < 311 ) %>% plotSpec()
betterData %>% filter( specNo==100, mass >1700, mass < 1702 ) %>% plotSpec()
# length(A)

# spec %>% filter(specNo==1) %>% select(mass,intensity) %>%
# bitonicClusters(1) %>% plot(type='l')


betterData %>% mutate( mass = round(mass/100)*100 ) %>% group_by( mass ) %>% summarise( intensity=sum(intensity) ) %>% plot



############################################################
# Enlighted bitonic clustering.


# betterData 	<- spec %>% split(.$specNo) %>% lapply(massDiffMass) %>% bind_rows(.id='specNo') 
# spec100 	<- betterData %>% filter(specNo==100, massDiff < 1. ) 
# chullIndices<- chull(spec100$mass, spec100$massDiff)
# pointsForFit<- spec100[chullIndices,]

polar <- 
	# betterData %>% filter( massDiff < .5, specNo==100 ) %>% mutate(
	betterData %>% filter(  specNo==100 ) %>% mutate(	
		sMd = sqrt(massDiff),
		Z 	= complex(real = mass, imaginary = sMd),
		arg = Arg(Z),
		mod = Mod(Z),
		tag = 'original'
	) 

polarHull <- 
	polar %>% bind_rows(.,
		data.frame(
			mass = 0,
			arg  = max(polar$arg)*10,
			tag  = 'dummy'
		)
	) %>% slice( chull(mass,arg) ) %>% 
	filter( tag != 'dummy' ) %>% # Adding a dummy point to make the convex hull as big as possible.
	mutate(mass2=mass^2)

polarHull %>% ggplot(aes(x=mass, y=arg))+geom_point()+theme_minimal()

# polarHull %>% ggplot(aes(x=mass,y=massDiff))+geom_point()+theme_minimal() 

LM <- polarHull %>% select(massDiff, mass, mass2) %>%  lm(formula = massDiff ~ mass + mass2)
predMassDiff <- predict( LM, data_clustered %>% mutate(mass2 = mass^2) %>% select(mass, mass2) ) 

plot(data_clustered %>% filter( massDiff < .5 ) %>% select(mass, massDiff),  pch='.')
lines(data_clustered$mass, predMassDiff, col='red')

pointsWithBiggerMassDiff <- data_clustered %>% mutate( res = massDiff-predMassDiff ) %>% filter( res > .01, massDiff<1., res < .1 ) 
pointsWithBiggerMassDiff %>% ggplot(aes(x=mass,y=massDiff,color=clustChange))+geom_point()
pointsWithBiggerMassDiff %>% ggplot(aes(x=mass,y=intensity,color=factor(cluster)))+geom_segment(aes(xend=mass,yend=0.0))+geom_point()+theme_minimal()+  theme(legend.position="none")
ggplotly()

pointsWithBiggerMassDiff %>% filter(mass>1382, mass < 1383)
testSpec %>% filter(mass>1375, mass < 1383) %>% ggplot(aes(x=mass,y=intensity,color=factor(cluster)))+geom_segment(aes(xend=mass,yend=0.0))+geom_point()+theme_minimal()+  theme(legend.position="none")
ggplotly()
# spec100 %>% filter( massDiff < .5 ) 


# polar$arg %>% summary

# polar %>% filter( arg < 10*median(arg) ) %>%

# polar %>% ggplot(aes(x=mass,y=arg))+geom_point()+theme_minimal() 
# # ggplot(aes(x=mass,y=arg))+geom_point()+theme_minimal() 


# LM 		<- pointsForFit %>% lm(formula = massDiff ~ I(mass) + I(mass^2), weights=1/(1+massDiff)) 
# LM 		<- pointsForFit[ LM$residuals < max(LM$residuals)/5,] %>% lm(formula = massDiff ~ mass + I(mass^2)) 
# coefs 	<- LM$coefficients * .8 # 1.25 is purely arbitrary.







# P <- spec100 %>% mutate(
# 	cluster = bitonicClustersEnlighted(spec100$mass, spec100$intensity, coefs ) %>% as.factor
# ) 

# P %>% filter(mass>1712,mass<1716) %>% 
# ggplot(aes(x=mass,y=intensity,color=cluster)) + geom_segment(aes(xend=mass,yend=0.0)) + theme_minimal()




# P %>% filter(mass>1712,mass<1716)


# # ggsave('longJohnPlot.pdf', plot=P, width=10000, height=4, limitsize=FALSE)
spec %>% filter( specNo == 1 ) %>% summary()
	
