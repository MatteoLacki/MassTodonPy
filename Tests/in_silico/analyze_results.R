spsm = suppressPackageStartupMessages
spsm(library("tidyverse"))

f = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/in_silico/results_Matteo/analyzed_100000.csv'
D = read_csv(f) %>% tbl_df
D = D %>% mutate( Q = factor(Q) ) 
levels(D$Q) = paste0('Q = ',levels(D$Q))  
  
P = D %>%
  ggplot(
  aes(x = corr, y=L1, color=eps)) +
  geom_point() +
  xlab('Correlation')+
  ylab('Normalized L1 Error (Max = 1)')+
  scale_colour_gradient(labels=scales::percent)+
  theme_minimal()+
  facet_grid(.~Q)
