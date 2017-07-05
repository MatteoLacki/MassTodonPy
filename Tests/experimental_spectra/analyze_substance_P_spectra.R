spsm = suppressPackageStartupMessages
spsm(library("jsonlite"))
spsm(library("tidyverse"))
len = length

analysis_results_path = '/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/experimental_spectra/analyzed_results.json'
R = read_json(analysis_results_path, simpifyVector=TRUE)

len(R)

parse_key = function(x) x %>% strsplit('_') %>% as.integer(x[-1]) %>% unlist

D = names(R) %>% strsplit('_') %>% 
    lapply(function(x) as.integer(x[-1])) %>% unlist %>% 
    matrix(ncol=3, byrow=TRUE) %>% 
    data.frame() %>% tbl_df %>% rename( Q=X1, WH=X2, WV=X3 )


lapply(
  R, 
  function(x) 
    
)
x = R[1]


parse_key(names(x[1]))

names(x)

lapply(x[[1]], function(y){
  y %>% unlist %>% data_frame(name=names(.),val=., ) 
})


R[[1]][[1]] %>% unlist %>% data_frame(name=names(.),val=.) 
R[[1]][[2]] %>% unlist %>% data_frame(name=names(.),val=.) 
R[[1]][[3]] %>% unlist %>% data_frame(name=names(.),val=.) 
