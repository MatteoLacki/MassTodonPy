library("jsonlite")
library("tidyverse")

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
# P(Y)
# data4plots$WH$parsedFiles
# P(data4plots$WH)

DATA = list()
i = 1
for(d_cnt in 1:length(data4plots)){
    X <- data4plots[[d_cnt]]$MassTodons
    for(exp_cnt in 1:length(X)){
        Y <- X[[exp_cnt]]
        DATA[[i]] = Y[c('precursor','maxProtonsNo')]
        DATA[[i]]$instrumental_setting = names(X)[exp_cnt]
        DATA[[i]]$mass_spectrum = getSpectrum( Y$spectrum )
        if( length(Y$cTerminus)>0 )
            DATA[[i]]$cTerminus = Y$cTerminus
        if( length(Y$nTerminus)>0 )
            DATA[[i]]$nTerminus = Y$nTerminus
        i = i+1
    }
}

storagePath = '/Users/matteo/Documents/MassTodon/Data/Storage/'
# pryr::object_size(toJSON(DATA,digits=10))
write_json(DATA, path='/Users/matteo/Documents/MassTodon/MassTodonPy/Tests/spectra.json', digits=10)
