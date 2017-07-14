id = 35
show_results = function(id){
  wrong = read.csv(paste0('/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/bootstrap/CSV_03_07_2017_mzPrec-065/', id,'.csv')) %>% 
          as_data_frame %>% 
          mutate(WH = factor(WH), WV = factor(WV) )
  res = list(
    'ETnoD' = wrong %>% select( real_or_bootstrap, contains('ETnoD')) ,
    'PTR'   = wrong %>% select( real_or_bootstrap, contains('PTR')) 
  )
  return(res)
}


lapply(res, View)
