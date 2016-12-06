library(dplyr);library(ggplot2);library(tidyr);library(Rcpp);library(plotly)
library(combinat)

# spec
load("/Users/matteo/Documents/Science/MassTodon/MassTodonPy/spectrumPreprocessing/FRL_220715_ubi_952_ETD_40ms_01.Rdata")

# bitonicLocalMatrix
sourceCpp('~/Documents/Science/MassTodon/MassTodonPy/spectrumPreprocessing/bitonicSequences.cpp')

# is the binning similar? 

spec <- spec %>% arrange(mass) 
# spec %>% group_by(specNo) %>% lapply(
# 	function(spectrum) bitonicLocal(spectrum,)
# )

# %>% mutate( clust = bitonicLocal(mass, intensity) ) 


spec100 <- spec %>% filter(specNo==100)
spec100 <- spec100 %>% mutate( clust=bitonicLocal(mass, intensity) ) %>% filter( intensity > 0 ) 

# spec100 %>% group_by( clust ) %>% summarise( n=n()) %>% plot(type='s')
spec100 %>% group_by( clust ) %>% summarise( n=n()) %>% filter(n==max(n))

plotSpec <- function(spectrum, k) spectrum %>% select(mass, intensity) %>% ggplot(aes(mass, intensity))+geom_line()+geom_point()
	
P <- spec100 %>% split(.$clust) %>% lapply(plotSpec)


normality <- function(mass,intensity){
	leftIdx <- which.min(mass)
	rightIdx<- which.max(mass)
	max(intensity) - intensity[leftIdx]  > abs(intensity[leftIdx]-intensity[rightIdx]) | 
	max(intensity) - intensity[rightIdx] > abs(intensity[leftIdx]-intensity[rightIdx])
}


# Finding the case where left is more similar to max than to right.
spec100 %>% group_by( clust ) %>% summarise( OK = normality(mass,intensity) ) %>% filter(!OK)


# spec100 %>% filter(clust==808)%>% select(mass, intensity) %>% ggplot(aes(mass, intensity))+geom_line()+geom_point()
# spec100 %>% filter(clust==809)%>% select(mass, intensity) %>% ggplot(aes(mass, intensity))+geom_line()+geom_point()
spec100 %>% filter(clust%in%c(808,809))%>% select(mass, intensity) %>% ggplot(aes(mass, intensity))+geom_line()+geom_point()
spec100 %>% filter(clust%in%c(808,809))%>% select(mass, intensity) %>% ggplot(aes(mass, intensity))+geom_line()+geom_point()


# spec100 %>% filter(clust==987)%>% select(mass, intensity) %>% ggplot(aes(mass, intensity))+geom_line()+geom_point()
# spec100 %>% filter(clust==988)%>% select(mass, intensity) %>% ggplot(aes(mass, intensity))+geom_line()+geom_point()
spec100 %>% filter(clust%in%c(987,988))%>% select(mass, intensity) %>% ggplot(aes(mass, intensity))+geom_line()+geom_point()
	# So, some peaks are not even approximately guassian-shaped.
	# This is all a big mess.

spec100 %>% filter(clust%in%c(984))%>% select(mass, intensity) %>% ggplot(aes(mass, intensity))+geom_line()+geom_point()
spec100 %>% filter(clust%in%c(983))%>% select(mass, intensity) %>% ggplot(aes(mass, intensity))+geom_line()+geom_point()



spec100$clust %>%  unique %>% max