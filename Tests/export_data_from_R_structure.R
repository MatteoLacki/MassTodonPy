suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("jsonlite"))

P = function(x, l=1) str(x,max.level=l)
getSpectrum = function(spec)
  spec %>% bind_rows %>% tbl_df %>% select(mz, intensity) %>% arrange(mz)

load('/Users/matteo/Dropbox/MassTodon/ProcessedData/SubstanceP/synapt_WH_WV_substanceP_010513.RData')

# P(data4plots)
# P(data4plots$WH)
# data4plots$WH$Data
# P(data4plots$WH$MassTodons)
# X <- data4plots$WH$MassTodons[[1]]
# P(X)
# Y <- X[[1]]
# Y$formulas
# Y$cTerminus
# length(Y$nTerminus)

getTerminalCnt <- function(terminus){
  if( length(terminus) > 0 )
      return(data.frame(terminus) %>% gather('Element','Count',1:ncol(.)))
  else
      return(NULL)
}

saveTerminalCnt <- function(terminus, path, i, type){
  term <- getTerminalCnt(terminus)
  if( !is.null(term) ){
      write_csv( term, path=paste0(path,i,'_', type,'_terminal.txt'), col_names=FALSE )
  }
}

# get other data and prepare a MetaInfo file with rows corresponding to i of a given spectrum.
storagePath = '/Users/matteo/Documents/MassTodon/Data/Storage/'

DATA <- list()

for(d_cnt in 1:length(data4plots)){
  X <- data4plots[[d_cnt]]$MassTodons
  for(exp_cnt in 1:length(X)){
    Y <- X[[exp_cnt]]
    mass_spectrum <- getSpectrum( Y$spectrum )
    write_csv( mass_spectrum, path=paste0(storagePath,i,'.txt'), col_names=FALSE )
  }
}
