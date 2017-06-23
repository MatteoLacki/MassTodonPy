library('tidyverse')

paths = c('data/probs_100000.csv',
          'data/probs_10000.csv')

d = Vectorize(function(a,b,A,B) sqrt((a-A)^2 + (b-B)^2)/sqrt(2) )

get_data = function(path){
  D = as_tibble(read.csv(path))
  D = D %>% mutate( rETnoD = ETnoD/(ETnoD+PTR), 
                    rPTR   = PTR/(ETnoD+PTR), 
                    int    = d(rETnoD, rPTR, ETnoD_int, PTR_int), 
                    base   = d(rETnoD, rPTR, ETnoD_base, PTR_base),
                    up_int = d(rETnoD, rPTR, ETnoD_up_int, PTR_up_int),
                    Q = factor(Q))
  levels(D$Q) = paste0('Q = ',levels(D$Q))
  return(D) }

D = lapply(paths, get_data) %>% bind_rows()
D$molsNo = factor(D$molsNo)
levels(D$molsNo) = c('10K IONS','100K IONS')

D %>% gather('analyzer', 'value', base, int, up_int ) %>%
  ggplot()+
  scale_y_continuous(labels=scales::percent) +
  geom_boxplot(aes(x=factor(Q),y=value,fill=analyzer)) +
  theme_minimal() +
  xlab('Initial Charge State')+
  ylab('Error') + 
  facet_grid(molsNo~.)
  
  