library(readMzXmlData)
library(dplyr)
library(scales)
library(ggplot2)
library(plotly)
library(tidyr)

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
	mzXmlFile 		= "FRL_220715_ubi_952_ETD_40ms_01.mzXML",
	removeMetaData 	= FALSE,
	verbose 		= TRUE
) 

spec <- lapply(
	Spectrum,
	function( info ) as.data.frame( info$spectrum )
) %>% bind_rows(.id='specNo') 	%>% 
filter(intensity > 10.) 		%>% tbl_df 

spec <- spec %>% mutate(specNo=as.integer(specNo))

metaDatas <- Spectrum %>% lapply(
	function(x)
		x$metaData[8:27] %>% data.frame
) %>% bind_rows %>% tbl_df

metaDatas$precursor.precursorCharge

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
	data.frame(
		mass 		= spectrum$mass[1:(N-1)],
		massDiff 	= sort(spectrum$mass) %>% diff 
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