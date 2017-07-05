source('get_data.R')

WH_150 = DD %>% filter(WH==150)

get_algo = Vectorize(function(string) strsplit(string,'_')[[1]][1]) 
check_if_precursor = Vectorize(function(string) tail(strsplit(string,'_')[[1]],1) == 'precursor') 

find_among_colnames(WH_150, 'basic')

